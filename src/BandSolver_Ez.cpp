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
 * iteration. Fortunately, this finite difference operator can be exactly
 * inverted in the Fourier domain.
 *  Alternatively, we can use a preconditioned Jacobi-Davidson method JDQZ
 * which requires a good preconditioner for (A - shift B), which the Fourier
 * domain inversion can be used as well.
 *
 * Details:
 *  Suppose we are in the compacted matrix indexing scheme, then we would like
 * to apply the FFT to each of the auxiliary fields in turn instead of all at
 * once to save memory. However, the inversion procedure is most straightforward
 * when all the FFT's are available. We note first that the submatrix containing
 * only P and V auxiliary variables is trivially inverted in the spatial domain,
 * while the electromagnetic field portion is easily inverted in the FFT domain.
 * The coupling between these two blocks is extremely low rank. Let us denote
 * the partitioned matrix system as:
 *     [ C  W^H ] [ u ] = [ g ]
 *     [ W  D   ] [ v ]   [ h ]
 * Left multiply to eliminate W^H:
 *     [ I   -W^H inv(D) ] [ C  W^H ] [ u ] = [ I   -W^H inv(D) ] [ g ]
 *     [ 0        I      ] [ W  D   ] [ v ]   [ 0        I      ] [ h ]
 * so
 *     [ C - W^H inv(D) W  0 ] [ u ] = [ g - W^H inv(D) h ]
 *     [         W         D ] [ v ]   [        h         ]
 * Also right multiply to eliminate W:
 *     [ C - W^H inv(D) W  0 ] [     I       0 ] [     I       0 ] [ u ] = [ g - W^H inv(D) h ]
 *     [         W         D ] [ -inv(D) W   I ] [  inv(D) W   I ] [ v ]   [        h         ]
 * giving
 *     [ C - W^H inv(D) W  0 ] [     I       0 ] [ u ] = [ g - W^H inv(D) h ]
 *     [         0         D ] [  inv(D) W   I ] [ v ]   [        h         ]
 * Precompute p = inv(D) h:
 *     [ C - W^H inv(D) W  0 ] [     I       0 ] [ u ] = [ g - W^H p ]
 *     [         0         I ] [  inv(D) W   I ] [ v ]   [    p      ]
 * Note that W^H inv(D) W has simple form:
 *   W^H inv(D) W = 
 * [    0     0     0     0 ] [   d1   -i wh1    0      0   ]^{-1} [ 0   i wp1   0 ]
 * [ -i wp1   0   i wp2   0 ] [  i wh1    q      0      0   ]      [ 0     0     0 ]
 * [    0     0     0     0 ] [    0      0      d2  -i wh2 ]      [ 0   i wp2   0 ]
 *                            [    0      0    i wh1    q   ]      [ 0     0     0 ]
 * = [ 0 0 0 ] where Q = q (wp1^2/det1 + wp2^2/det2)
 *   [ 0 Q 0 ]   where det = d q - wh^2
 *   [ 0 0 0 ]     Note that we encounter problems when d q = wh^2.
 * This represents a diagonal modification of the E field in real space, which
 * is simple to apply with an FFT. Therefore, we proceed to compute
 *   u = inv(C - W^H inv(D) W) (g - W^H p)
 * Therefore, we finally have an algorithm; compute in order:
 *   v = inv(D) h
 *   u = inv(C - W^H inv(D) W) (g - W^H v)
 *   v += inv(D) W u
 * When d q = wh^2, then Q becomes infinite, corresponding to choosing a shift
 * exactly equal to a resonance frequency. When this happens, the fields are
 * null inside the material. We can forcibly zero out the fields inside the
 * material, or we can try to avoid it.
 *
 * Roadmap to implementing general meshes:
 *  Build adjacency tables which give the coefficients of the curl operators
 *   Matrix construction loops over number of unknowns per grid cell
 *  Build corresponding Hodge star ratios (single constants for uniform mesh)
 *  Build analogous FFT operator, and invert.
 *
 */

typedef std::complex<double> complex_t;

#define DIVH_OFF 3
#define EX_OFF   (undefined)
#define EY_OFF   (undefined)
#define EZ_OFF   1
#define HX_OFF   2
#define HY_OFF   0
#define HZ_OFF   (undefined)

#define IDX(i,j) (res[1]*(i)+(j))
#define UNIDX(q,i,j) do{ \
	(i) = (q) / res[1]; \
	(j) = (q) % res[1]; \
	}while(0)

SPB::BandSolver_Ez::BandSolver_Ez(double Lr[4]):BandSolver(Lattice2(Lr)),L(Lr),
	N(0),B(NULL),
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
				ind2cell[2*cell2ind[q]+1] = q;
				
				// get materials of this cell (for now, simple pointwise check)
				matind[q] = 0;
				int tag, num_poles;
				double p[2] = {
					L.Lr[0]*fi + L.Lr[2]*fj,
					L.Lr[1]*fi + L.Lr[3]*fj
				};
				if(!shapeset.QueryPt(p, &tag)){
					tag = -1;
				}
				if(-1 == tag){
					num_poles = 0;
				}else{
					num_poles = material[tag].poles.size();
					used_mat.insert(tag);
				}
				npoles[q] = num_poles;
				if(tag >= 0){
					matind[q] = tag+1; // assume it is within 4 bit limit
				}
				
			}
		}
		
		int next_pole_offset = 0;
		for(std::set<size_t>::const_iterator i = used_mat.begin(); i != used_mat.end(); ++i){
			mat_to_pole_offset[*i] = next_pole_offset;
			next_pole_offset += material[*i].poles.size();
		}
		
		// Now perform the indexing, need to set cell2ind and ind2cell[2*p+0];
		int next_index = 0;
		for(int p = 0; p < Ngrid; ++p){
			int q = ind2cell[2*p+1];
			int i,j; UNIDX(q,i,j);
			ind2cell[2*p+0] = next_index;
			cell2ind[q] = next_index;
			
			// We have Ez, Hx, Hy, divH, plus V,P pairs for each pole
			next_index += 4 + 2*npoles[q];
		}
		N = next_index;
		if(AtGamma){
			constraint_off = N;
			n_constraint = 4 + 2*next_pole_offset;
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
					
					const int curmat = matind[q];
					complex_t eps_z(1.);
					if(curmat > 0){
						eps_z = material[curmat-1].eps_inf.value[8];
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
				row = row0 + 4;
				int matbits = matind[q];
				while(matbits & 0xF){
					const int m = (matbits & 0xF)-1;
					complex_t eps_z(1.);
					if(m > 0){
						eps_z = material[m-1].eps_inf.value[8];
					}
					for(int pi = 0; pi < material[m].poles.size(); ++pi){
						const LorentzPole &pole = material[m].poles[pi];
						if(0 == pole.omega_0){
							ASET(row,row, complex_t(0.));
							ASET(row+1,row, complex_t(0., 1.));
							ASET(row,row+1, complex_t(0.,-1.));
							ASET(row+1,row+1, -target);
							BSET(row,row, 0.);
						}else{
							ASET(row,row, -target);
							ASET(row+1,row, complex_t(0., pole.omega_0) * eps_z);
							ASET(row,row+1, complex_t(0.,-pole.omega_0) * eps_z);
							ASET(row+1,row+1, -target);
							BSET(row,row, 1.);
						}
						ASET(row0+EZ_OFF, row+1, complex_t(0.,-pole.omega_p) * eps_z);
						ASET(row+1, row0+EZ_OFF, complex_t(0., pole.omega_p) * eps_z);
						BSET(row+1,row+1, 1.);
						
						if(constraint_off){
							int constraint_row = constraint_off + 4 + mat_to_pole_offset[m] + 2*pi;
							ASET(constraint_row, row+0, complex_t(0, -1));
							ASET(row+0, constraint_row, complex_t(0,  1));
							
							constraint_row++;
							ASET(constraint_row, row+1, complex_t(0, -1));
							ASET(row+1, constraint_row, complex_t(0,  1));
						}
						
						row += 2;
					}
					matbits >>= 4;
				}
			}
		}
		for(int i = 0; i < n_constraint; ++i){
			BSET(constraint_off+i,constraint_off+i, 0.);
			ASET(constraint_off+i,constraint_off+i, 10.);
			if(0 == i%2){
				ASET(constraint_off+i+1,constraint_off+i, complex_t(0,0));
			}else{
				ASET(constraint_off+i-1,constraint_off+i, complex_t(0,0));
			}
		}
	}
	B = new sparse_t(N,N, Bmap);
	Atmp = new sparse_t(N,N, Amap);
	if(0){
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
	for(int j = 0; j < 4*N; ++j){
		double *d = (double*)ldl_data.D;
		printf("%g\n", d[j]);
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
