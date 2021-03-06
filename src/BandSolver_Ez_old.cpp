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
	N(0),B(NULL),Atmp(NULL),
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
	
	// delete A
    if(valid_A_numeric){
		free(ldl_data.Lp);
		free(ldl_data.parent);
		free(ldl_data.Li);
		free(ldl_data.Lx);
		free(ldl_data.D);
	}
	if(NULL != B){ delete B; }
	if(NULL != Atmp){ delete Atmp; }
}

size_t SPB::BandSolver_Ez::GetSize() const{
	return N;
}
void SPB::BandSolver_Ez::Aop(const std::complex<double> *x, std::complex<double> *y) const{
	RNP::Sparse::MultMV<'N'>(*Atmp, x, y);
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
	
	if(valid_A_numeric){
		free(ldl_data.Lp);
		free(ldl_data.parent);
		free(ldl_data.Li);
		free(ldl_data.Lx);
		free(ldl_data.D);
	}
	if(NULL != B){
		delete B;
	}
	if(NULL != Atmp){
		delete Atmp;
	}
	return 0;
}

void Acol1(int col, int *nnz, int *row_ind, void *data){
	RNP::Sparse::TCCSMatrix<complex_t> *A = (RNP::Sparse::TCCSMatrix<complex_t>*)data;
	const int k1 = A->colptr[2*col+0];
	const int k2 = A->colptr[2*col+1];
	const int nnz1 = k2-k1;
	const int nnz2 = A->colptr[2*col+2] - k2;
	int p1 = 0;
	int p2 = 0;
//printf("col = %d\n", col);
	*nnz = 0;
	for(; p1 < nnz1 && p2 < nnz2; ){
		int i1 = A->rowind[k1+p1] / 2;
		int i2 = A->rowind[k2+p2] / 2;
//printf(" i1,2 = %d,%d\n", i1, i2);
		if(i1 < i2){
			if(i1 >= col){ break; }
		}else{
			if(i2 >= col){ break; }
		}
		
		if(i1 < i2){
			*row_ind = i1;
			row_ind++;
			(*nnz)++;
			p1++;
		}else if(i1 > i2){
			*row_ind = i2;
			row_ind++;
			(*nnz)++;
			p2++;
		}else{
			*row_ind = i1;
			row_ind++;
			(*nnz)++;
			p1++;
			p2++;
		}
	}
}
void Acol2(int col, int *nnz, int *row_ind, double *row_val, void *data){
	RNP::Sparse::TCCSMatrix<complex_t> *A = (RNP::Sparse::TCCSMatrix<complex_t>*)data;
	const int k1 = A->colptr[2*col+0];
	const int k2 = A->colptr[2*col+1];
	const int nnz1 = k2-k1;
	const int nnz2 = A->colptr[2*col+2] - k2;
	int p1 = 0;
	int p2 = 0;
//printf("col = %d\n", col);
	
	*nnz = 0;
	for(; p1 < nnz1 && p2 < nnz2;){
		int i1raw = A->rowind[k1+p1];
		int i1 = i1raw / 2;
		int i1off = i1raw % 2;
		int i2raw = A->rowind[k2+p2];
		int i2 = i2raw / 2;
		int i2off = i2raw % 2;
		if(i1 < i2){
			if(i1 > col){ break; }
		}else{
			if(i2 > col){ break; }
		}
		
		int p1inc = 0;
		int p2inc = 0;
		
		for(int j = 0; j < 8; ++j){ row_val[j] = 0; }
		
		if(i1 < i2){
			*row_ind = i1;
			row_val[2*(0+i1off)+0] = A->values[k1+p1].real();
			row_val[2*(0+i1off)+1] = A->values[k1+p1].imag();
//printf("A(%d,%d) = ", i1, col); C2PRINT(row_val);
			row_ind++;
			row_val += 8;
			(*nnz)++;
			p1++;
		}else if(i1 > i2){
			*row_ind = i2;
			row_val[2*(2+i2off)+0] = A->values[k2+p2].real();
			row_val[2*(2+i2off)+1] = A->values[k2+p2].imag();
//printf("A(%d,%d) = ", i2, col); C2PRINT(row_val);
			row_ind++;
			row_val += 8;
			(*nnz)++;
			p2++;
		}else{
			*row_ind = i1;
			if(i1 == col){ // diag is full
				row_val[2*(0)+0] = A->values[k1+p1].real();
				row_val[2*(0)+1] = A->values[k1+p1].imag();
				row_val[2*(1)+0] = A->values[k1+p1+1].real();
				row_val[2*(1)+1] = A->values[k1+p1+1].imag();
				row_val[2*(2)+0] = A->values[k2+p2].real();
				row_val[2*(2)+1] = A->values[k2+p2].imag();
				row_val[2*(3)+0] = A->values[k2+p2+1].real();
				row_val[2*(3)+1] = A->values[k2+p2+1].imag();
//printf("A(%d,%d) = ", i1, col); C2PRINT(row_val);
			}else{
				row_val[2*(0+i1off)+0] = A->values[k1+p1].real();
				row_val[2*(0+i1off)+1] = A->values[k1+p1].imag();
				row_val[2*(2+i2off)+0] = A->values[k2+p2].real();
				row_val[2*(2+i2off)+1] = A->values[k2+p2].imag();
//printf("A(%d,%d) = ", i1, col); C2PRINT(row_val);
			}
			row_ind++;
			row_val += 8;
			(*nnz)++;
			p1++;
			p2++;
			if(i1 == col){ break; }
		}
		p1 += p1inc;
		p2 += p2inc;
	}fflush(stdout);
}


int SPB::BandSolver_Ez::MakeASymbolic(){
	const size_t Ngrid = res[0] * res[1];
	
	const double klen = hypot(last_k[0], last_k[1]);
	const bool AtGamma = (klen < std::numeric_limits<double>::epsilon() * L.CharacteristicKLength());
	
	if(NULL != ind2cell){ free(ind2cell); }
	if(NULL != cell2ind){ free(cell2ind); }
	if(NULL != matind){ free(matind); }
	if(NULL != npoles){ free(npoles); }
	ind2cell = (int*)malloc(sizeof(int) * 2*Ngrid);
	cell2ind = (int*)malloc(sizeof(int) * Ngrid);
	matind = (int*)malloc(sizeof(int) * Ngrid);
	npoles = (int*)malloc(sizeof(int) * Ngrid);

	// Prepare the indexing
	std::map<size_t,size_t> mat_to_pole_offset;
	std::set<size_t> used_mat;
	int constraint_off = 0, n_constraint = 0;
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
						//if(0 != curmat.poles[pi].omega_0){
							num_poles++;
						//}
					}
					used_mat.insert(tag);
				}
				npoles[q] = num_poles;
				matind[q] = tag+2; // assume it is within 4 bit limit
			}
		}
		
		int next_pole_offset = 0;
		for(std::set<size_t>::const_iterator i = used_mat.begin(); i != used_mat.end(); ++i){
			mat_to_pole_offset[*i] = next_pole_offset;
			const Material& curmat = material[*i];
			for(int pi = 0; pi < curmat.poles.size(); ++pi){
				next_pole_offset++;
				//if(0 != curmat.poles[pi].omega_0){
					next_pole_offset++;
				//}
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
		if(AtGamma){
			constraint_off = N;
			n_constraint = NUM_EH + 2*NUM_E*next_pole_offset;
			N += n_constraint;
		}
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
					
					if(constraint_off){
						ASET(constraint_off+HY_OFF, col, complex_t(0, 1));
						ASET(col, constraint_off+HY_OFF, complex_t(0,-1));
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
					
					if(constraint_off){
						ASET(constraint_off+DIVH_OFF, col, complex_t(0, 1));
						ASET(col, constraint_off+DIVH_OFF, complex_t(0,-1));
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
					
					if(constraint_off){
						ASET(constraint_off+EZ_OFF, col, complex_t(0, 1));
						ASET(col, constraint_off+EZ_OFF, complex_t(0,-1));
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
					
					if(constraint_off){
						ASET(constraint_off+HX_OFF, col, complex_t(0, 1));
						ASET(col, constraint_off+HX_OFF, complex_t(0,-1));
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
							int V_off = 1;
							
							if(0 != pole.omega_0){
								ASET(row,row, -target*eps_z);
								ASET(row+V_off,row, complex_t(0., pole.omega_0) * eps_z);
								ASET(row,row+V_off, complex_t(0.,-pole.omega_0) * eps_z);
								BSET(row,row, eps_z);
//std::cout << "P!=0, V row = " << row+V_off << std::endl;
							}else{
								ASET(row,row, 0.);
								ASET(row+V_off,row, complex_t(0., 1));
								ASET(row,row+V_off, complex_t(0.,-1));
								BSET(row,row, 0.);
//std::cout << "P=0, V row = " << row+V_off << std::endl;
							}
							ASET(row0+EZ_OFF, row+V_off, complex_t(0.,-pole.omega_p) * eps_z);
							ASET(row+V_off, row0+EZ_OFF, complex_t(0., pole.omega_p) * eps_z);
							ASET(row+V_off,row+V_off, -target*eps_z);
							BSET(row+V_off,row+V_off, eps_z);
						
							if(constraint_off){
								int constraint_row = constraint_off + NUM_EH + mat_to_pole_offset[m] + 2*pi;
								if(0 != pole.omega_0){
									ASET(constraint_row, row+0, complex_t(0, -1));
									ASET(row+0, constraint_row, complex_t(0,  1));
								}
								constraint_row++;
								ASET(constraint_row, row+V_off, complex_t(0, -1));
								ASET(row+V_off, constraint_row, complex_t(0,  1));
							}
						
							row += 2;
						}
					}
					matbits >>= 4;
				}
			}
		}
		for(int i = 0; i < n_constraint; ++i){
			BSET(constraint_off+i,constraint_off+i, 0.);
			ASET(constraint_off+i,constraint_off+i, 0.);
			if(0 == i%2){
				ASET(constraint_off+i+1,constraint_off+i, complex_t(0,0));
			}else{
				ASET(constraint_off+i-1,constraint_off+i, complex_t(0,0));
			}
		}
	}
	B = new sparse_t(N,N, Bmap);
	Atmp = new sparse_t(N,N, Amap);
	if(1){
		std::cout << "A="; RNP::Sparse::PrintSparseMatrix(*Atmp) << ";" << std::endl;
		std::cout << "B="; RNP::Sparse::PrintSparseMatrix(*B) << ";" << std::endl;
		exit(0);
	}
	
	ldl_data.Lp = (int*)malloc(sizeof(int) * (N+1));
	ldl_data.parent = (int*)malloc(sizeof(int) * N);
	ldl_data.D  = (double*)malloc(sizeof(double) * 4*N);
	
	int *Lnz = (int*)malloc(sizeof(int) * N);
	int *Flag = (int*)malloc(sizeof(int) * N);
	int *Pattern = (int*)malloc(sizeof(int) * N);
	double *Y = (double*)malloc(sizeof(double) * 4*N);
	int *rowind = (int*)malloc(sizeof(int) * N);
	double *rowval = (double*)malloc(sizeof(double) * 8*N);
	LDL_symbolic(N/2, N, &Acol1, ldl_data.Lp, ldl_data.parent, Lnz, Flag, rowind, (void*)Atmp);
	const int lnz = ldl_data.Lp[N/2];
	ldl_data.Li = (int*)malloc(sizeof(int) * lnz);
	ldl_data.Lx = (double*)malloc(sizeof(double) * 8*lnz);
	LDL_numeric(N/2, N, &Acol2, ldl_data.Lp, ldl_data.parent, Lnz, ldl_data.Li, ldl_data.Lx, ldl_data.D, Y, Pattern, Flag, rowind, rowval, (void*)Atmp);
	
	// It would appear that the matrix remains singular (well, the D matrix)
	// even with all the extra constraints. This only seems to happen at
	// zero shift, so we just have to avoid a target frequency of zero at
	// the Gamma point.
	
	/*
	if(AtGamma){ // last part needs zeroing
		double *d = (double*)ldl_data.D;
		for(int j = 0; j < n_constraint; ++j){
			int k = 4*N-4*n_constraint+4*j;
			d[k+0] = 0;
			d[k+1] = 0;
			d[k+2] = 0;
			d[k+3] = 0;
		}
	}*/
	/*
	for(int j = 0; j < N/2; j+=8){
		double *d = (double*)ldl_data.D;
		C2PRINT(&d[j]);
	}
	*/
	free(rowval);
	free(rowind);
	free(Y);
	free(Pattern);
	free(Flag);
	free(Lnz);
	
	valid_A_numeric = true;
	
	return 0;
}

void SPB::BandSolver_Ez::ShiftInv(const complex_t &shift, const complex_t *x, complex_t *y) const{
	int info;
	RNP::TBLAS::Copy(N, x,1, y,1);
	
	LDL_lsolve(N/2, (double*)y, ldl_data.Lp, ldl_data.Li, ldl_data.Lx);
	LDL_dsolve(N/2, (double*)y, ldl_data.D);
	LDL_ltsolve(N/2, (double*)y, ldl_data.Lp, ldl_data.Li, ldl_data.Lx);
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

size_t SPB::BandSolver_Ez::GetProblemSize() const{
	return N;
}
