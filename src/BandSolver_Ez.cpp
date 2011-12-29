#define _USE_MATH_DEFINES
#define RNP_OUTPUT_MATHEMATICA
//#define RNP_OUTPUT_MATLAB
#define RNP_SPARSE_USE_IO

#include "SPB.hpp"
#include <cstdlib>
#include <limits>
#include <ctime>
#include <IO.h>
#include <Sparse.h>
#include <cstring>
#include <cstdio>
#include <fstream>
#include "util.h"
#include <LinearSolve.h>
#include <set>
extern "C"{
#include "NestedDissection.h"
#include "ldl2.h"
}

/* Problem setup:
 *  In Ez polarization, we have fields Hx, Hy, and Ez.
 * We discretize the fundamental parallelogram with a primal-dual
 * mesh. For the square lattice, this corresponds to the flattened
 * Yee-cell:
 *
 *     +-------------+
 *     |             |
 *     |     divH    |     ^ y
 *  Hx ->     +      |     |
 *     |             |     +---> x
 *     |      ^      |
 *     +------|------+
 *    Ez      Hy
 *
 * The Ez are located on the primal grid vertices, while the H
 * components are located on the dual grid edges. The div(H)
 * constraint is enforced on dual vertices.
 *  For a non-orthogonal cell, we have two cases:
 *
 *         +-----------+
 *        / \         /
 *       /   \_  +   /         _
 *  Hv <-    |\ Hw  /          /| v
 *     /   +   \   /          /        u.v > 0
 *    /    ^    \ /          +---> u
 *   +-----|-----+
 *  Ez    Hu
 *
 * The case with u.v < 0 is similar:
 *
 *     +-----------+
 *      \         / \
 *       \   +  _/   \        __ v
 *    Hv <-     /|HW  \      |\
 *         \   /   +   \       \       u.v < 0
 *          \ /    ^    \       +---> u
 *           +-----|-----+
 *           Ez    Hu
 *
 * We solve the eigenvalue problem (for example with 2 resonances):
 *    [    0       0     i grad     0      0      0      0     ] [dvH] (div H)
 *    [    0       0     i curl  -i wp1    0    i wp2    0     ] [ E ]
 *    [ -i grad -i curl     0       0      0      0      0     ] [ H ]
 *    [    0     i wp1      0    -i G1  -i w01    0      0     ] [ V1]
 *    [    0       0        0     i w01    0      0      0     ] [ P1]
 *    [    0     i wp2      0       0      0   -i G2  -i w02   ] [ V2]
 *    [    0       0        0       0      0    i w01    0     ] [ P2]
 *  = diag( 0, eps, mu, 1, 1, 1, 1) (dvH,E,H,V1,P1,V2,P2)
 * where wp is the plasmon frequency \omega_p times \epsilon_\infty
 *   and eps is \epsilon_\infty
 *   and eta is 1/\epsilon_\infty
 *   and w0 is \omega_0 times \epsilon_\infty
 *   and G is \Gamma times \epsilon_\infty
 *  The ordering of matrix entries is the same as above, except each of
 * the variables H (Hx and Hy), E, and dvH are blocks of size Ngrid = Nx*Ny.
 * The remaining auxiliary fields are indexed in blocks for each grid cell.
 * The index offsets from the end of dvH are stored in the array `ind`.
 *
 *  Writing this out as a matrix system, we get
 *     (A - lambda B) x = 0
 * where A is full rank, and B is rank deficient and approximately diagonal.
 * When A and B are Hermitian, the eigenvalues come in real pairs +/-lambda.
 * with eigenvectors that are complex conjugates of each other.
 *  We have a number of ways to solve the eigensystem. Since we are typically
 * interested in eigenvalues in the interior of the spectrum, we can use a
 * spectral transformation:
 *   inv(A - shift B) B x = mu x      mu = 1/(lambda - shift)
 * This requires solving the linear system (A - shift B) x = b for each outer
 * iteration.
 *
 * Roadmap to implementing general meshes:
 *  Build adjacency tables which give the coefficients of the curl operators
 *   Matrix construction loops over number of unknowns per grid cell
 *  Build corresponding Hodge star ratios (single constants for uniform mesh)
 *
 */

typedef std::complex<double> complex_t;

#define NDIM 2
#define NUM_EH   4
#define NUM_E    1
#define DIVE_OFF (undefined)
#define EX_OFF   (undefined)
#define HZ_OFF   (undefined)
#define EY_OFF   (undefined)
#define HY_OFF   0
#define EZ_OFF   1
#define HX_OFF   2
#define DIVH_OFF 3

/*
#define NLATVECS 3
// NLATVECS gives number of rows of the following table, which always has 3 columns
// Edir below references the row of this table.
static const int lattice_tab[] = {
	1,0,0,
	0,1,0,
	0,0,1
}
static const int field_tab[] = {
	// field, Bval, Nterms, Edir
	// -------------------------
	// field = which field is on the LHS
	// Bval = 0 means B is 0
	//        1 means B is 1
	//        2 means B is epsilon
	// Nterms is number of finite differences to other EH fields (number of coupling to other EH fields divided by 2)
	// Edir = -1 if not an E-field
	//      = 0-n for lattice vector direction (possibly more than 3 for skew lattices)
	
	// Hy
	HY_OFF  , 1, 2, -1, // Nterms = 2 for Ez, and 2 for divH
	// Ez
	EZ_OFF  , 2, 2,  2, // Nterms = 2 each for Hx,Hy
	// Hx
	HX_OFF  , 1, 2, -1, // Nterms = 2 for Ez, and 2 for divH
	// divH
	DIVH_OFF, 0, 2, -1
};
static const int diff_tab[] = {
	// diff dim, fb, sign, field
	// diff dim = direction in which to take finite difference (x=0 or y=1)
	// fb = forward=0 or backward=1 difference
	// sign = sign of the difference
	// field = field to difference (use one of the _OFF values)
	
	// Hy
	0, 0, 0, EZ_OFF,   // idr[0] * (Ez[i+1,j,k] - Ez[i,j,k])
	1, 1, 0, DIVH_OFF, // idr[1] * (divH[i,j,k] - divH[i,j-1,k])
	// Ez
	1, 1, 1, HX_OFF, // -idr[1] * (Hx[i,j,k] - Hx[i,j-1,k])
	0, 1, 0, HY_OFF, //  idr[0] * (Hy[i,j,k] - Hy[i-1,j,k])
	// Hx
	1, 0, 1, EZ_OFF,  // -idr[1] * (Ez[i,j+1,k] - Ez[i,j,k])
	0, 1, 0, DIVH_OFF //  idr[0] * (divH[i,j,k] - divH[i-1,j,k])
	// divH
	0, 0, 0, HX_OFF, // idr[0] * (Hx[i+1,j,k] - Hx[i,j,k])
	1 ,0, 0, HY_OFF  // idr[1] * (Hy[i,j+1,k] - Hy[i,j,k])
};

void table_driven(){
	for(ii[0] = 0; ii[0] < res[0]; ++ii[0]){
		for(ii[1] = 0; ii[1] < res[1]; ++ii[1]){
			const int q = IDX(ii);
			const int row0 = cell2ind[q];
			int row = 0, col;
			
			// Get Epsilon tensor for E-H field coupling
			complex_t epsE[9];
			// Get Epsilon tensor for E-V field coupling
			complex_t epsV[9];
			
			int diff_tab_offset = 0;
			for(int field = 0; field < NUM_EH; ++field){
				row = row0 + field_tab[4*field+0];
				switch(field_tab[4*field+1]){
				case 0: BSET(col,col, 0); break;
				case 1: BSET(col,col, 1); break;
				default:
					// set block of B for E-field directions
					break;
				}
				
				// set the diagonal contribution
				//ASET(col,col, -target * Bvals);
				
				static const double sign[2] = {1., -1.};
				
				for(int diff = 0; diff < field_tab[4*field+2]; ++diff, ++diff_tab_offset){
					const int *curdiff = &diff_tab[4*diff_tab_offset];
					const int dim = curdiff[0];
					const complex_t coeff = complex_t(0,sign[curdiff[2]] * idr[dim]);
					if(curdiff[1]){ // backward diff
						int jj[NDIM]; for(int _i=0;_i<NDIM;++_i){jj[_i]=ii[_i];}
						if(0 == ii[dim]){
							col = row0 + curdiff[3];
							ASET(row,col, coeff);
							jj[dim] = res[dim]-1;
							col = cell2ind[IDX(jj)] + curdiff[3];
							ASET(row,col, -coeff/Bloch[dim]);
						}else{
							jj[dim]--;
							row = cell2ind[IDX(jj)] + curdiff[3];
							ASET(row,col, -coeff);
							col = row0 + curdiff[3];
							ASET(row,col, coeff);
						}
					}else{
						if(ii[dim] + 1 == res[dim]){
						}else{
						}
					}
				}
				
				if(constraint_off){
					// add constraints
				}
			}
			
			// For each E-field, add corresponding P and V fields
		}
	}
}*/

#define IDX(i,j) (res[1]*(i)+(j))
#define UNIDX(q,i,j) do{ \
	(i) = (q) / res[1]; \
	(j) = (q) % res[1]; \
	}while(0)

SPB::BandSolver_Ez::BandSolver_Ez(double Lr[4]):BandSolver(Lattice2(Lr)),L(Lr),
	N(0),B(NULL),A(NULL),ldl(),
	valid_A_numeric(false),
	last_shift(0)
{
	cell2ind = NULL;
	ind2cell = NULL;
	matind = NULL;
	npoles = NULL;
}

SPB::BandSolver_Ez::~BandSolver_Ez(){
	
	free(cell2ind);
	free(ind2cell);
	free(matind);
	free(npoles);
	if(NULL != B){ delete B; }
	if(NULL != A){ delete A; }
}

size_t SPB::BandSolver_Ez::GetProblemSize() const{
	return N;
}
size_t SPB::BandSolver_Ez::GetSize() const{
	return N;
}
void SPB::BandSolver_Ez::Aop(const std::complex<double> *x, std::complex<double> *y) const{
	RNP::Sparse::MultMV<'N'>(*A, x, y);
}
void SPB::BandSolver_Ez::Bop(const std::complex<double> *x, std::complex<double> *y) const{
	RNP::Sparse::MultMV<'N'>(*B, x, y);
}
void SPB::BandSolver_Ez::PrintField(const std::complex<double> *y, const char *filename) const{
	const int Ngrid = res[0]*res[1];
	std::ostream *out = &std::cout;
	if(NULL != filename){
		out = new std::ofstream(filename);
	}
	/*
	for(int i = 0; i < res[0]; ++i){
		for(int j = 0; j < res[1]; ++j){
			(*out)
				<< i << "\t" << j
				 << "\t" << y[HX_OFF + IDX(i,j)].real()
				 << "\t" << y[HX_OFF + IDX(i,j)].imag()
				 << "\t" << y[HY_OFF + IDX(i,j)].real()
				 << "\t" << y[HY_OFF + IDX(i,j)].imag()
				 << "\t" << y[EZ_OFF + IDX(i,j)].real()
				 << "\t" << y[EZ_OFF + IDX(i,j)].imag()
				 << "\t" << y[DIVH_OFF+IDX(i,j)].real()
				 << "\t" << y[DIVH_OFF+IDX(i,j)].imag()
				 << "\n";
		}
		(*out) << "\n";
	}*/
	if(NULL != filename){
		delete out;
	}
}

void SPB::BandSolver_Ez::Precond(const std::complex<double> &alpha, const std::complex<double> &beta, const std::complex<double> *x, std::complex<double> *y) const{
	ShiftInv(alpha/beta, x, y);
	return;
}
void SPB::BandSolver_Ez::Orth(std::complex<double> *y) const{
}
void SPB::BandSolver_Ez::Randvec(std::complex<double> *y) const{
}
	
int SPB::BandSolver_Ez::SolveK(const double *k){
	SPB_VERB(1, "Solving k-point (%.14g, %.14g)\n", k[0], k[1]);
	ClearSolution();
	
	last_k[0] = k[0];
	last_k[1] = k[1];
	UpdateA(k);
	IRASolve(N);
	
	last_k[0] = k[0];
	last_k[1] = k[1];
	return 0;
}

int SPB::BandSolver_Ez::OutputEpsilon(int *res, const char *filename, const char *format) const{
	if(NULL == res){ return -1; }
	if(NULL != format && 0 == strcmp(format, "gnuplot")){
		FILE *fp = stdout;
		if(NULL != filename){
			fp = fopen(filename, "wb");
		}
		if(2 == dim){
			double p[2];
			for(int i = 0; i < res[0]; ++i){
				double fi = ((double)i) / (double)res[0] - 0.5;
				for(int j = 0; j < res[1]; ++j){
					double fj = ((double)j) / (double)res[1] - 0.5;
					p[0] = L.Lr[0]*fi + L.Lr[2]*fj;
					p[1] = L.Lr[1]*fi + L.Lr[3]*fj;
					int tag;
					if(!shapeset.QueryPt(p, &tag)){ tag = -1; }
					fprintf(fp, "%d\t%d\t%d\n", i, j, tag);
				}
				fprintf(fp, "\n");
			}
		}else{
		}
		if(NULL != filename){
			fclose(fp);
		}
	}
	return 0;
}









	
int SPB::BandSolver_Ez::InvalidateByStructure(){
	
	free(cell2ind); cell2ind = NULL;
	free(ind2cell); ind2cell = NULL;
	free(matind); matind = NULL;
	free(npoles); npoles = NULL;
	
	if(NULL != B){
		delete B;
	}
	return 0;
}


int SPB::BandSolver_Ez::GetMaxBlockSize() const{
	return 2;
}
int SPB::BandSolver_Ez::GetMaxNNZPerRow() const{
	return N;
}
int SPB::BandSolver_Ez::BeginBlockSymbolic() const{
	assembly_data.i = 0;
	assembly_data.j = 0;
	assembly_data.k = 0;
}
int SPB::BandSolver_Ez::BeginBlockNumeric() const{
	assembly_data.i = 0;
	assembly_data.j = 0;
	assembly_data.k = 0;
}
int SPB::BandSolver_Ez::GetNextBlockSymbolic(int *rowptr, int *colind) const{
	if(assembly_data.i >= N){ return 0; }
	rowptr[0] = A->colptr[assembly_data.i];
	rowptr[1] = A->colptr[assembly_data.i+1];
	rowptr[2] = A->colptr[assembly_data.i+2];
	for(int r = 0; r < 2; ++r){
		for(int p = rowptr[r]; p < rowptr[r+1]; ++p){
			*colind = A->rowind[p];
			colind++;
		}
	}
	rowptr[1] -= rowptr[0];
	rowptr[2] -= rowptr[0];
	rowptr[0] = 0;
	assembly_data.i += 2;
	return 2;
}
int SPB::BandSolver_Ez::GetNextBlockNumeric(int *rowptr, int *colind, complex_t *value) const{
	if(assembly_data.i >= N){ return 0; }
	rowptr[0] = A->colptr[assembly_data.i];
	rowptr[1] = A->colptr[assembly_data.i+1];
	rowptr[2] = A->colptr[assembly_data.i+2];
	for(int r = 0; r < 2; ++r){
		for(int p = rowptr[r]; p < rowptr[r+1]; ++p){
			*colind = A->rowind[p];
			colind++;
			*value = A->values[p];
			value++;
		}
	}
	rowptr[1] -= rowptr[0];
	rowptr[2] -= rowptr[0];
	rowptr[0] = 0;
	assembly_data.i += 2;
	return 2;
}

int SPB::BandSolver_Ez::MakeASymbolic(){
	const size_t Ngrid = res[0] * res[1];
	
	const double klen = hypot(last_k[0], last_k[1]);
	
	if(NULL != ind2cell){ free(ind2cell); }
	if(NULL != cell2ind){ free(cell2ind); }
	if(NULL != matind){ free(matind); }
	if(NULL != npoles){ free(npoles); }
	ind2cell = (int*)malloc(sizeof(int) * 2*Ngrid);
	cell2ind = (int*)malloc(sizeof(int) * Ngrid);
	matind = (int*)malloc(sizeof(int) * Ngrid);
	npoles = (int*)malloc(sizeof(int) * Ngrid);

	// Prepare the indexing
	std::set<size_t> used_mat;
	{
		{ // use cell2ind to temporarily hold the forward permutation
			int a[2] = {res[1], 1};
			int p[2] = {1,1};
			NestedDissectionCartesian1v(2, res, a, p, cell2ind);
		}
		for(int i = 0; i < res[0]; ++i){
			const double fi = ((double)i/(double)res[0]) - 0.5;
			for(int j = 0; j < res[1]; ++j){
				const double fj = ((double)j/(double)res[1]) - 0.5;
				const int q = IDX(i,j);
				//ind2cell[2*cell2ind[q]+1] = q;
				ind2cell[2*q+1] = q;
				
				// get materials of this cell (for now, simple pointwise check)
				matind[q] = 0;
				int tag, num_poles = 0;
				double p[2] = {
					L.Lr[0]*fi + L.Lr[2]*fj,
					L.Lr[1]*fi + L.Lr[3]*fj
				};
				if(!shapeset.QueryPt(p, &tag)){
					tag = -1;
				}
				if(tag >= 0){
					const Material& curmat = material[tag];
					for(int pi = 0; pi < curmat.poles.size(); ++pi){
						num_poles++;
						num_poles++;
					}
					used_mat.insert(tag);
				}
				npoles[q] = num_poles;
				matind[q] = tag+2; // assume it is within 4 bit limit
			}
		}
		
		int next_pole_offset = 0;
		for(std::set<size_t>::const_iterator i = used_mat.begin(); i != used_mat.end(); ++i){
			const Material& curmat = material[*i];
			for(int pi = 0; pi < curmat.poles.size(); ++pi){
				next_pole_offset++;
				next_pole_offset++;
			}
		}
		
		// Now perform the indexing, need to set cell2ind and ind2cell[2*p+0];
		int next_index = 0;
		for(int p = 0; p < Ngrid; ++p){
			int q = ind2cell[2*p+1];
			int i,j; UNIDX(q,i,j);
			ind2cell[2*p+0] = next_index;
			cell2ind[q] = next_index;
			
			// We have Ez, Hx, Hy, divH, plus V,P pairs for each pole
			next_index += NUM_EH + NUM_E*npoles[q];
		}
		N = next_index;
	}
	
	sparse_t::entry_map_t Amap;
	sparse_t::entry_map_t Bmap;
	
	//complex_t *Adata = (complex_t*)doublecomplexMalloc(Annz);
	//int *rowind = intMalloc(Annz);
	//int *colptr = intMalloc((N+1));
	
	{
		//size_t Aind = 0;
		const double Lrl[2] = {
			hypot(L.Lr[0], L.Lr[1]),
			hypot(L.Lr[2], L.Lr[3])
		};
		const double idr[2] = {
			(double)res[0] / Lrl[0],
			(double)res[1] / Lrl[1]
		};
		
		double use_k[2] = { last_k[0], last_k[1] };
				
		const complex_t Bloch[2] = {
			complex_t(cos(use_k[0]*2*M_PI), sin(use_k[0]*2*M_PI)),
			complex_t(cos(use_k[1]*2*M_PI), sin(use_k[1]*2*M_PI))
		};
		
/*
		Adata[Aind] = (COEFF); \
		rowind[Aind] = (ROW); \
		Aind++; \
*/
#define ASET(ROW,COL,COEFF) do{ \
		Amap[sparse_t::index_t((ROW),(COL))] =  (COEFF); \
	}while(0)
//Amap[sparse_t::index_t((COL),(ROW))] =  std::conj(COEFF);
#define ASETCOL(COL,IND) /*colptr[(COL)] = (IND)*/
#define BSET(ROW,COL,COEFF) Bmap[sparse_t::index_t((ROW),(COL))] = (COEFF)

		size_t row;
		complex_t coeff;

		for(int i = 0; i < res[0]; ++i){
			for(int j = 0; j < res[1]; ++j){
				const int q = IDX(i,j);
				const int row0 = cell2ind[q];
				int row = 0, col;
				
				{ // Set the Hy column
					col = row0 + HY_OFF;
					BSET(col,col, 1);
					ASET(col,col, -target);
					
					// Ez = complex_t(0,-idr[1]) * (Hx[i,j,k] - Hx[i,j-1,k])
					//    + complex_t(0, idr[0]) * (Hy[i,j,k] - Hy[i-1,j,k])
					coeff = complex_t(0,idr[0]);
					if(0 == i){
						row = row0 + EZ_OFF;
						ASET(row,col, coeff);
						row = cell2ind[IDX(res[0]-1,j)] + EZ_OFF;
						ASET(row,col, -coeff/Bloch[0]);
					}else{
						row = cell2ind[IDX(i-1,j)] + EZ_OFF;
						ASET(row,col, -coeff);
						row = row0 + EZ_OFF;
						ASET(row,col, coeff);
					}
					
					// divH = idr[0] * (Hx[i+1,j,k] - Hx[i,j,k])
					//      + idr[1] * (Hy[i,j+1,k] - Hy[i,j,k]) <--
					coeff = complex_t(0,idr[1]);
					if(j+1 == res[0]){
						row = cell2ind[IDX(i,0)] + DIVH_OFF;
						ASET(row,col, coeff*Bloch[1]);
						row = row0 + DIVH_OFF;
						ASET(row,col, -coeff);
					}else{
						row = row0 + DIVH_OFF;
						ASET(row,col, -coeff);
						row = cell2ind[IDX(i,j+1)] + DIVH_OFF;
						ASET(row,col, coeff);
					}
				}
				
				{ // Set the divH column
					int col = row0 + DIVH_OFF;
					ASETCOL(col,Aind);
					BSET(col,col, 0);
					ASET(col,col, complex_t(0.));

					// H += i grad divH
					// Hx += idr[0] * (divH[i,j,k] - divH[i-1,j,k])
					// Hy += idr[1] * (divH[i,j,k] - divH[i,j-1,k])
					
					coeff = complex_t(0,idr[1]);
					if(0 == j){
						row = row0 + HY_OFF;
						ASET(row,col, coeff);
						row = cell2ind[IDX(i,res[1]-1)] + HY_OFF;
						ASET(row,col, -coeff/Bloch[1]);
					}else{
						row = cell2ind[IDX(i,j-1)] + HY_OFF;
						ASET(row,col, -coeff);
						row = row0 + HY_OFF;
						ASET(row,col, coeff);
					}
					
					coeff = complex_t(0,idr[0]);
					if(0 == i){
						row = row0 + HX_OFF;
						ASET(row,col, coeff);
						row = cell2ind[IDX(res[0]-1,j)] + HX_OFF;
						ASET(row,col, -coeff/Bloch[0]);
					}else{
						row = cell2ind[IDX(i-1,j)] + HX_OFF;
						ASET(row,col, -coeff);
						row = row0 + HX_OFF;
						ASET(row,col, coeff);
					}
				}
				
				{ // Set the Ez column
					col = row0 + EZ_OFF;
					ASETCOL(col,Aind);
					
					const int curmat = matind[q]-2;
					complex_t eps_z(1.);
					if(curmat >= 0){
						eps_z = material[curmat].eps_inf.value[8];
					}
					BSET(col,col, eps_z);

					ASET(col,col, -target*eps_z);
					
					// Hx = complex_t(0,-idr[1]) * (Ez[i,j+1,k] - Ez[i,j,k])
					coeff = complex_t(0,idr[1]);
					if(j+1 == res[1]){
						row = cell2ind[IDX(i,0)] + HX_OFF;
						ASET(row,col, -coeff*Bloch[1]);
						row = row0 + HX_OFF;
						ASET(row,col, coeff);
					}else{
						row = row0 + HX_OFF;
						ASET(row,col, coeff);
						row = cell2ind[IDX(i,j+1)] + HX_OFF;
						ASET(row,col, -coeff);
					}
					
					// Hy = complex_t(0, idr[0]) * (Ez[i+1,j,k] - Ez[i,j,k])
					coeff = complex_t(0,idr[0]);
					if(i+1 == res[0]){
						row = cell2ind[IDX(0,j)] + HY_OFF;
						ASET(row,col, coeff*Bloch[0]);
						row = row0 + HY_OFF;
						ASET(row,col, -coeff);
					}else{
						row = row0 + HY_OFF;
						ASET(row,col, -coeff);
						row = cell2ind[IDX(i+1,j)] + HY_OFF;
						ASET(row,col, coeff);
					}
				}
				
				{ // Set the Hx column
					int col = row0 + HX_OFF;
					BSET(col,col, 1);
					
					// divH = idr[0] * (Hx[i+1,j,k] - Hx[i,j,k]) <--
					//      + idr[1] * (Hy[i,j+1,k] - Hy[i,j,k])
					coeff = complex_t(0,idr[0]);
					if(i+1 == res[0]){
						row = cell2ind[IDX(0,j)] + DIVH_OFF;
						ASET(row,col, coeff*Bloch[0]);
						row = row0 + DIVH_OFF;
						ASET(row,col, -coeff);
					}else{
						row = row0 + DIVH_OFF;
						ASET(row,col, -coeff);
						row = cell2ind[IDX(i+1,j)] + DIVH_OFF;
						ASET(row,col, coeff);
					}
					
					// Ez = complex_t(0,-idr[1]) * (Hx[i,j,k] - Hx[i,j-1,k])
					//    + complex_t(0, idr[0]) * (Hy[i,j,k] - Hy[i-1,j,k])
					coeff = complex_t(0,idr[1]);
					if(0 == j){
						row = row0 + EZ_OFF;
						ASET(row,col, -coeff);
						row = cell2ind[IDX(i,res[1]-1)] + EZ_OFF;
						ASET(row,col, coeff/Bloch[1]);
					}else{
						row = cell2ind[IDX(i,j-1)] + EZ_OFF;
						ASET(row,col, coeff);
						row = row0 + EZ_OFF;
						ASET(row,col, -coeff);
					}

					ASET(col,col, -target);
				}
				
				// Deal with all the poles
				// oh god that sounds racist
				row = row0 + NUM_EH;
				int matbits = matind[q];
				while(matbits & 0xF){
					int m = (matbits & 0xF)-1;
					complex_t eps_z(1.);
					if(m > 0){
//std::cout << i << "\t" << j << "\t" << m << std::endl;
						m--;
						eps_z = material[m].eps_inf.value[8];
//std::cout << i << "\t" << j << "\t" << eps_z << std::endl;
						const Material& curmat = material[m];
						for(int pi = 0; pi < curmat.poles.size(); ++pi){
							const LorentzPole &pole = material[m].poles[pi];
//std::cout << "\t" << pi << "\t" << pole.omega_0 << "\t" << pole.omega_p << std::endl;

						
							ASET(row,row, -target*eps_z);
							ASET(row+1,row, complex_t(0., pole.omega_0) * eps_z);
							ASET(row,row+1, complex_t(0.,-pole.omega_0) * eps_z);
							BSET(row,row, eps_z);
//std::cout << "P!=0, V row = " << row+V_off << std::endl;
							ASET(row0+EZ_OFF, row+1, complex_t(0.,-pole.omega_p) * eps_z);
							ASET(row+1, row0+EZ_OFF, complex_t(0., pole.omega_p) * eps_z);
							ASET(row+1,row+1, -target*eps_z);
							BSET(row+1,row+1, eps_z);
						
							row += 2;
						}
					}
					matbits >>= 4;
				}
			}
		}
	}
	B = new sparse_t(N,N, Bmap);
	A = new sparse_t(N,N, Amap);
	if(0){
		std::cout << "A="; RNP::Sparse::PrintSparseMatrix(*A) << ";" << std::endl;
		std::cout << "B="; RNP::Sparse::PrintSparseMatrix(*B) << ";" << std::endl;
		exit(0);
	}
	
	ldl.Factorize(*this);
	
	valid_A_numeric = true;
	
	return 0;
}

void SPB::BandSolver_Ez::ShiftInv(const complex_t &shift, const complex_t *x, complex_t *y) const{
	int info;
	RNP::TBLAS::Copy(N, x,1, y,1);
	
	ldl.Solve(y);
	/*
	complex_t *res = new complex_t[N];
	double *Y = (double*)res;
	for(int i = 0; i < N; i++){
		res[i] = -x[i];
	}
	RNP::Sparse::MultMV<'N'>(*Atmp, y, res, complex_t(1.), complex_t(1.));
	
	//for(int j = 0; j < N/2; j++){
	//	printf("Y[%d] = {%.14g+i%.14g %.14g+i%.14g}\n", j, Y[4*j+0], Y[4*j+1], Y[4*j+2], Y[4*j+3]);
	//}
	// rnorm = norm (y, inf)
	double rnorm = 0;
	for(int i = 0; i < 2*N; i++){
		double r = (Y[i] > 0) ? (Y[i]) : (-Y[i]);
		rnorm = (r > rnorm) ? (r) : (rnorm);
	}
	printf ("relative maxnorm of residual: %g\n", rnorm);
	
	delete [] res;
	*/
	//std::cout << "info  = " << info << std::endl;
}

int SPB::BandSolver_Ez::MakeANumeric(){
	return 0;
}
int SPB::BandSolver_Ez::UpdateA(const double k[2]){
	//temp hacky
	InvalidateByStructure();
	if(NULL == cell2ind){
		MakeASymbolic();
	}
	if(!valid_A_numeric){
		MakeANumeric();
	}
	if(last_k[0] == k[0] && last_k[1] == k[1]){
		return 0;
	}
	// do update
	return 0;
}

int SPB::BandSolver_Ez::GetN() const{
	return N;
}
