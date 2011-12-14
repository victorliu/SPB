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

#define DIVH_OFF (0*Ngrid)
#define EX_OFF   (undefined)
#define EY_OFF   (undefined)
#define EZ_OFF   (1*Ngrid)
#define HX_OFF   (2*Ngrid)
#define HY_OFF   (3*Ngrid)
#define HZ_OFF   (undefined)

#define IDX(i,j) (res[1]*(i)+(j))

SPB::BandSolver_Ez::BandSolver_Ez(double Lr[4]):BandSolver(Lattice2(Lr)),L(Lr),
	N(0),ind(NULL),B(NULL),
	valid_A_numeric(false),
	last_shift(0)
{
    set_default_options(&superlu_data.options);
    StatInit(&superlu_data.stat);
}

SPB::BandSolver_Ez::~BandSolver_Ez(){
	free(ind);
	// delete A
    if(valid_A_numeric){
		if(NULL != superlu_data.perm_r){ SUPERLU_FREE(superlu_data.perm_r); }
		if(NULL != superlu_data.perm_c){ SUPERLU_FREE(superlu_data.perm_c); }
		Destroy_CompCol_Matrix(&superlu_data.A);
		Destroy_SuperNode_Matrix(&superlu_data.L);
		Destroy_CompCol_Matrix(&superlu_data.U);
	}
	if(NULL != B){ delete B; }
	if(NULL != Atmp){ delete Atmp; }
    StatFree(&superlu_data.stat);
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










void
zgsfact(superlu_options_t *options, SuperMatrix *A, int *perm_c, int *perm_r,
      SuperMatrix *L, SuperMatrix *U,
      SuperLUStat_t *stat, int *info ){
    SuperMatrix AC; /* Matrix postmultiplied by Pc */
    int      lwork = 0, *etree, i;
    
    /* Set default values for some parameters */
    int      panel_size;     /* panel size */
    int      relax;          /* no of columns in a relaxed snodes */
    int      permc_spec;

    *info = 0;

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = NATURAL:  natural ordering 
     *   permc_spec = MMD_AT_PLUS_A: minimum degree on structure of A'+A
     *   permc_spec = MMD_ATA:  minimum degree on structure of A'*A
     *   permc_spec = COLAMD:   approximate minimum degree column ordering
     *   permc_spec = MY_PERMC: the ordering already supplied in perm_c[]
     */
    permc_spec = COLAMD;
    if ( permc_spec != MY_PERMC && options->Fact == DOFACT )
      get_perm_c(permc_spec, A, perm_c);

    etree = intMalloc(A->ncol);

    sp_preorder(options, A, perm_c, etree, &AC);

    panel_size = sp_ienv(1);
    relax = sp_ienv(2);

    /*printf("Factor PA = LU ... relax %d\tw %d\tmaxsuper %d\trowblk %d\n", 
	  relax, panel_size, sp_ienv(3), sp_ienv(4));*/
    /* Compute the LU factorization of A. */
    
    zgstrf(options, &AC, relax, panel_size, etree,
            NULL, lwork, perm_c, perm_r, L, U, stat, info);

    SUPERLU_FREE(etree);
    Destroy_CompCol_Permuted(&AC);
}

	
int SPB::BandSolver_Ez::InvalidateByStructure(){
	if(NULL != ind){ free(ind); ind = NULL; }
	if(valid_A_numeric){
		if(NULL != superlu_data.perm_r){ SUPERLU_FREE(superlu_data.perm_r); }
		if(NULL != superlu_data.perm_c){ SUPERLU_FREE(superlu_data.perm_c); }
		Destroy_SuperNode_Matrix(&superlu_data.L);
		Destroy_CompCol_Matrix(&superlu_data.U);
		Destroy_CompCol_Matrix(&superlu_data.A);
	}
	if(NULL != B){
		delete B;
	}
	return 0;
}
int SPB::BandSolver_Ez::MakeASymbolic(){
	const size_t Ngrid = res[0] * res[1];
	
	const double klen = hypot(last_k[0], last_k[1]);
	const bool AtGamma = (klen < std::numeric_limits<double>::epsilon() * L.CharacteristicKLength());
	
	size_t Annz = 0;
	
	size_t extra_constraints = 0;
	size_t EH_constraints = 0;
	
//#define ONLY_LOWER
	// First block column: divH, couples to Hx and Hy twice each + diag
	Annz += 5*Ngrid;
	// Second block column: Ez, couples to Hx and Hy twice each + diag
	Annz += 5*Ngrid;
	// Each block of Hx and Hy contributes a diagonal
	Annz += 2*Ngrid;
#ifndef ONLY_LOWER
	// Each block of Hx and Hy couples to divH and Ez twice each
	Annz += 2*4*Ngrid;
#endif
	// Each row of a V contributes 2 nonzeros + 1 possible diaonal
	// Each row of a P contributes 1 nonzeros + 1 possible diaonal
	if(NULL != ind){ free(ind); }
	ind = (int*)malloc(sizeof(int) * 2*Ngrid);

	// Prepare the indexing
	size_t extra_constraint_start = 0;
	size_t EH_constraint_start = 0;
	std::map<size_t,size_t> used_mat_poles;
	std::map<size_t,size_t> mat_counts;
	{
		size_t next_index = 0;
		for(int i = 0; i < res[0]; ++i){
			const double fi = ((double)i/(double)res[0]) - 0.5;
			for(int j = 0; j < res[1]; ++j){
				const double fj = ((double)j/(double)res[1]) - 0.5;
				ind[2*IDX(i,j)+0] = 4*Ngrid+next_index;
				
				// get material of this cell (simple pointwise check)
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
//std::cout << i << "\t" << j << "\t" << impl->eps_z_fft[IDX(i,j)] << "\t" << num_poles << std::endl;
				}
				ind[2*IDX(i,j)+1] = tag;
				// update next index
				size_t zero_poles = 0;
				for(size_t p = 0; p < num_poles; ++p){
					if(0 == material[tag].poles[p].omega_0){
						zero_poles++;
					}
				}
				next_index += 2*(num_poles-zero_poles);
			
				Annz += (num_poles-zero_poles) * 4;
#ifndef ONLY_LOWER
				Annz += (num_poles-zero_poles) * 2;
#endif
				mat_counts[tag]++;
				if(used_mat_poles[tag] == 0){
					used_mat_poles[tag] = num_poles-zero_poles;
				}
			}
		}
		if(AtGamma){
			size_t last_offset = extra_constraints; // should be 0
			for(std::map<size_t,size_t>::iterator i = used_mat_poles.begin(); i != used_mat_poles.end(); ++i){
				size_t n_poles = i->second;
				i->second = last_offset;
				extra_constraints += n_poles;
				Annz += n_poles * (mat_counts[i->first]+1);
	#ifndef ONLY_LOWER
				Annz += n_poles * mat_counts[i->first];
	#endif
				last_offset = n_poles;
			}
		}
		extra_constraint_start = 4*Ngrid + next_index;
		
		EH_constraint_start = extra_constraint_start+extra_constraints;
		/*
		// for divH
		EH_constraints += 1;
			Annz += 1*(Ngrid + 1);
#ifndef ONLY_LOWER
			Annz += 1*Ngrid;
#endif*/

		if(AtGamma){
			EH_constraints += 3;
			Annz += 3*(Ngrid + 1);
#ifndef ONLY_LOWER
			Annz += 3*Ngrid;
#endif
		}
		
		N = 4*Ngrid + next_index + extra_constraints + EH_constraints;
	}
//std::cout << "extra_constraints = " << extra_constraints << std::endl;
	
	sparse_t::entry_map_t Amap;
	sparse_t::entry_map_t Bmap;
	
	complex_t *Adata = (complex_t*)doublecomplexMalloc(Annz);
	int *rowind = intMalloc(Annz);
	int *colptr = intMalloc((N+1));
	
	{
		size_t Aind = 0;
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
		
#define ASET(ROW,COL,COEFF) do{ \
		Adata[Aind] = (COEFF); \
		rowind[Aind] = (ROW); \
		Aind++; \
		Amap[sparse_t::index_t((ROW),(COL))] =  (COEFF); \
	}while(0)
//Amap[sparse_t::index_t((COL),(ROW))] =  std::conj(COEFF);
#define ASETCOL(COL,IND) colptr[(COL)] = (IND)
#define BSET(ROW,COL,COEFF) Bmap[sparse_t::index_t((ROW),(COL))] = (COEFF)

		size_t row;
		complex_t coeff;
		
		for(int i = 0; i < res[0]; ++i){
			for(int j = 0; j < res[1]; ++j){
				const size_t col = DIVH_OFF + IDX(i,j);
				ASETCOL(col,Aind);
				BSET(col,col, 0);
				ASET(col,col, complex_t(0.));

				// H += i grad divH
				// Hx += idr[0] * (divH[i,j,k] - divH[i-1,j,k])
				// Hy += idr[1] * (divH[i,j,k] - divH[i,j-1,k])
				coeff = complex_t(0,idr[0]);
				if(0 == i){
					row = HX_OFF + IDX(i,j);
					ASET(row,col, coeff);
					row = HX_OFF + IDX(res[0]-1,j);
					ASET(row,col, -coeff/Bloch[0]);
				}else{
					row = HX_OFF + IDX(i-1,j);
					ASET(row,col, -coeff);
					row = HX_OFF + IDX(i,j);
					ASET(row,col, coeff);
				}
				
				coeff = complex_t(0,idr[1]);
				if(0 == j){
					row = HY_OFF + IDX(i,j);
					ASET(row,col, coeff);
					row = HY_OFF + IDX(i,res[1]-1);
					ASET(row,col, -coeff/Bloch[1]);
				}else{
					row = HY_OFF + IDX(i,j-1);
					ASET(row,col, -coeff);
					row = HY_OFF + IDX(i,j);
					ASET(row,col, coeff);
				}
				if(EH_constraints){
					ASET(EH_constraint_start+0,col, complex_t(1.));
				}
			}
		}
		for(int i = 0; i < res[0]; ++i){
			for(int j = 0; j < res[1]; ++j){
				const size_t col = EZ_OFF + IDX(i,j);
				ASETCOL(col,Aind);
				
				const int curmat = ind[2*IDX(i,j)+1];
				complex_t eps_z(1.);
				if(curmat >= 0){
					eps_z = material[curmat].eps_inf.value[8];
				}
				BSET(col,col, eps_z);
				ASET(col,col, -target*eps_z);
				
				// Hx = complex_t(0,-idr[1]) * (Ez[i,j+1,k] - Ez[i,j,k])
				coeff = complex_t(0,idr[1]);
				if(j+1 == res[1]){
					row = HX_OFF + IDX(i,0);
					ASET(row,col, -coeff*Bloch[1]);
					row = HX_OFF + IDX(i,j);
					ASET(row,col, coeff);
				}else{
					row = HX_OFF + IDX(i,j);
					ASET(row,col, coeff);
					row = HX_OFF + IDX(i,j+1);
					ASET(row,col, -coeff);
				}
				
				// Hy = complex_t(0, idr[0]) * (Ez[i+1,j,k] - Ez[i,j,k])
				coeff = complex_t(0,idr[0]);
				if(i+1 == res[0]){
					row = HY_OFF + IDX(0,j);
					ASET(row,col, coeff*Bloch[0]);
					row = HY_OFF + IDX(i,j);
					ASET(row,col, -coeff);
				}else{
					row = HY_OFF + IDX(i,j);
					ASET(row,col, -coeff);
					row = HY_OFF + IDX(i+1,j);
					ASET(row,col, coeff);
				}
				
				
				if(curmat >= 0){
					const Material &m = material[curmat];
					const size_t np = m.poles.size();
					size_t po = 0; // actual pole offset (ignoring poles @ 0)
					const int row0 = ind[2*IDX(i,j)+0];
					for(size_t p = 0; p < np; ++p){
						if(0 == m.poles[p].omega_0){ continue; }
						row = row0 + 2*po + 0; // V_p
				
						coeff = complex_t(0, m.poles[p].omega_p) * eps_z;
						ASET(row,col, coeff);
						
						po++;
					}
				}
			}
		}
		
		for(int i = 0; i < res[0]; ++i){
			for(int j = 0; j < res[1]; ++j){
				const size_t col = HX_OFF + IDX(i,j);
				BSET(col,col, 1);
				ASETCOL(col,Aind);
#ifndef ONLY_LOWER
				
				// divH = idr[0] * (Hx[i+1,j,k] - Hx[i,j,k]) <--
				//      + idr[1] * (Hy[i,j+1,k] - Hy[i,j,k])
				coeff = complex_t(0,idr[0]);
				if(i+1 == res[0]){
					row = DIVH_OFF + IDX(0,j);
					ASET(row,col, coeff*Bloch[0]);
					row = DIVH_OFF + IDX(i,j);
					ASET(row,col, -coeff);
				}else{
					row = DIVH_OFF + IDX(i,j);
					ASET(row,col, -coeff);
					row = DIVH_OFF + IDX(i+1,j);
					ASET(row,col, coeff);
				}
				
				// Ez = complex_t(0,-idr[1]) * (Hx[i,j,k] - Hx[i,j-1,k])
				//    + complex_t(0, idr[0]) * (Hy[i,j,k] - Hy[i-1,j,k])
				coeff = complex_t(0,idr[1]);
				if(0 == j){
					row = EZ_OFF + IDX(i,j);
					ASET(row,col, -coeff);
					row = EZ_OFF + IDX(i,res[1]-1);
					ASET(row,col, coeff/Bloch[1]);
				}else{
					row = EZ_OFF + IDX(i,j-1);
					ASET(row,col, coeff);
					row = EZ_OFF + IDX(i,j);
					ASET(row,col, -coeff);
				}
#endif
				ASET(col,col, -target);
				if(EH_constraints > 1){
					ASET(EH_constraint_start+1,col, complex_t(1.));
				}
			}
		}
		for(int i = 0; i < res[0]; ++i){
			for(int j = 0; j < res[1]; ++j){
				const size_t col = HY_OFF + IDX(i,j);
				BSET(col,col, 1);
				ASETCOL(col,Aind);
#ifndef ONLY_LOWER
				
				// divH = idr[0] * (Hx[i+1,j,k] - Hx[i,j,k])
				//      + idr[1] * (Hy[i,j+1,k] - Hy[i,j,k]) <--
				coeff = complex_t(0,idr[1]);
				if(j+1 == res[0]){
					row = DIVH_OFF + IDX(i,0);
					ASET(row,col, coeff*Bloch[1]);
					row = DIVH_OFF + IDX(i,j);
					ASET(row,col, -coeff);
				}else{
					row = DIVH_OFF + IDX(i,j);
					ASET(row,col, -coeff);
					row = DIVH_OFF + IDX(i,j+1);
					ASET(row,col, coeff);
				}
				
				// Ez = complex_t(0,-idr[1]) * (Hx[i,j,k] - Hx[i,j-1,k])
				//    + complex_t(0, idr[0]) * (Hy[i,j,k] - Hy[i-1,j,k])
				coeff = complex_t(0,idr[0]);
				if(0 == i){
					row = EZ_OFF + IDX(i,j);
					ASET(row,col, coeff);
					row = EZ_OFF + IDX(res[0]-1,j);
					ASET(row,col, -coeff/Bloch[0]);
				}else{
					row = EZ_OFF + IDX(i-1,j);
					ASET(row,col, -coeff);
					row = EZ_OFF + IDX(i,j);
					ASET(row,col, coeff);
				}
#endif
				ASET(col,col, -target);
				if(EH_constraints > 1){
					ASET(EH_constraint_start+2,col, complex_t(1.));
				}
			}
		}
		for(int i = 0; i < res[0]; ++i){
			for(int j = 0; j < res[1]; ++j){
				const int col0 = ind[2*IDX(i,j)+0];
				size_t col;
				const int curmat = ind[2*IDX(i,j)+1];
				complex_t eps_z(1.);
				if(curmat >= 0){
					eps_z = material[curmat].eps_inf.value[8];
				}
//std::cerr << "starting Aind = " << Aind;
				if(curmat >= 0){
					const Material &m = material[curmat];
					const size_t np = m.poles.size();
					size_t po = 0; // actual pole offset (ignoring poles @ 0)
					for(size_t p = 0; p < np; ++p){
						if(0 == m.poles[p].omega_0){ continue; }
						col = col0 + 2*po + 0; // V_p
						ASETCOL(col,Aind);
#ifndef ONLY_LOWER
						row = EZ_OFF + IDX(i,j); // E
						coeff = complex_t(0, m.poles[p].omega_p) * eps_z;
						ASET(row,col, -coeff);
#endif
						coeff = complex_t(0,-m.poles[p].Gamma) * eps_z;
						ASET(col,col, coeff - target);
						BSET(col,col, 1);
						
						coeff = complex_t(0, m.poles[p].omega_0) * eps_z;
						row = col0 + 2*po + 1; // P
						ASET(row,col, coeff);
						
						col = col0 + 2*po + 1; // P
						row = col0 + 2*po + 0; // V

						ASETCOL(col,Aind);
						BSET(col,col, 1);
#ifndef ONLY_LOWER
						ASET(row,col, -coeff);
#endif
						ASET(col,col, -target);
						if(extra_constraints){
//std::cerr << "used_mat_poles[curmat] = " << used_mat_poles[curmat] << ", po = " << po << std::endl;
							ASET(extra_constraint_start+used_mat_poles[curmat]+po,col, complex_t(1.));
						}
						po++;
					}
				}
//std::cerr << ", ending Aind = " << Aind << std::endl;
			}
		}
		if(extra_constraints){
			size_t col = extra_constraint_start;
			for(std::map<size_t,size_t>::const_iterator m = used_mat_poles.begin(); m != used_mat_poles.end(); ++m){
				const Material &mat = material[m->first];
				const size_t np = mat.poles.size();
				size_t po = 0; // actual pole offset (ignoring poles @ 0)
				for(size_t p = 0; p < np; ++p){
					if(0 == mat.poles[p].omega_0){ continue; }
						
					ASETCOL(col,Aind);
					BSET(col,col, 0);
#ifndef ONLY_LOWER
					for(int i = 0; i < res[0]; ++i){
						for(int j = 0; j < res[1]; ++j){
							if(m->first == ind[2*IDX(i,j)+1]){
								row = ind[2*IDX(i,j)+0] + 2*po + 1; // P
						
								ASET(row,col, complex_t(1.));
//std::cerr << "Setting row = " << row << ", col = " << col << std::endl;
							}
						}
					}
#endif
					ASET(col,col, complex_t(0.));
					po++;
					col++;
				}
			}
		}
		//std::cerr << "EH_constraints = " << EH_constraints << " extra_constraints = " << extra_constraints << std::endl;
		//std::cerr << "EH_constraint_start = " << EH_constraint_start << " extra_constraint_start = " << extra_constraint_start << std::endl;
		if(EH_constraints > 0){
			size_t col = EH_constraint_start+0;
			ASETCOL(col,Aind);
			BSET(col,col, 0);
#ifndef ONLY_LOWER
			for(int i = 0; i < res[0]; ++i){
				for(int j = 0; j < res[1]; ++j){
					const size_t row = DIVH_OFF + IDX(i,j);
					ASET(row,col, complex_t(1.));
				}
			}
#endif
			ASET(col,col, complex_t(0.));
		}
		if(EH_constraints > 1){
			size_t col = EH_constraint_start+1;
			ASETCOL(col,Aind);
			BSET(col,col, 0);
#ifndef ONLY_LOWER
			for(int i = 0; i < res[0]; ++i){
				for(int j = 0; j < res[1]; ++j){
					const size_t row = HX_OFF + IDX(i,j);
					ASET(row,col, complex_t(1.));
				}
			}
#endif
			ASET(col,col, complex_t(0.));
			
			col = EH_constraint_start+2;
			ASETCOL(col,Aind);
			BSET(col,col, 0);
#ifndef ONLY_LOWER
			for(int i = 0; i < res[0]; ++i){
				for(int j = 0; j < res[1]; ++j){
					const size_t row = HY_OFF + IDX(i,j);
					ASET(row,col, complex_t(1.));
				}
			}
#endif
			ASET(col,col, complex_t(0.));
		}
		
		ASETCOL(N,Aind);
		//std::cerr << "Aind = " << Aind << ", Annz = " << Annz << std::endl;
	}
	B = new sparse_t(N,N, Bmap);
	Atmp = new sparse_t(N,N, Amap);
	if(0){
		std::cout << "A="; RNP::Sparse::PrintSparseMatrix(*Atmp) << ";" << std::endl;
		std::cout << "B="; RNP::Sparse::PrintSparseMatrix(*B) << ";" << std::endl;
		exit(0);
	}
	zCreate_CompCol_Matrix(&superlu_data.A, N, N, Annz,
		(doublecomplex*)Adata, rowind, colptr, SLU_NC, SLU_Z,
#ifdef ONLY_LOWER
		SLU_SYL
#else
		SLU_GE
#endif
		);
	
	if(0){
		NCformat *ncf = (NCformat*)superlu_data.A.Store;
		for(int i = 0; i <= Atmp->m; ++i){
			std::cout << Atmp->colptr[i] << "\t" << ncf->colptr[i];
			if(Atmp->colptr[i] != ncf->colptr[i]){
				std::cout << "\t*";
			}
			std::cout << std::endl;
		}
		for(int i = 0; i < ncf->nnz; ++i){
			std::cout << Atmp->rowind[i] << "\t" << ncf->rowind[i];
			if(Atmp->rowind[i] != ncf->rowind[i]){
				std::cout << "\t*";
			}
			std::cout << std::endl;
		}
		exit(0);
	}
	
	int info;
	superlu_data.perm_c = intMalloc(N);
	superlu_data.perm_r = intMalloc(N);
	zgsfact(&superlu_data.options, &superlu_data.A,
		superlu_data.perm_c, superlu_data.perm_r,
		&superlu_data.L, &superlu_data.U,
		&superlu_data.stat, &info);
	//std::cout << "info = " << info << std::endl;
	valid_A_numeric = true;
	return 0;
}

void SPB::BandSolver_Ez::ShiftInv(const complex_t &shift, const complex_t *x, complex_t *y) const{
	int info;
	SuperMatrix B;
	RNP::TBLAS::Copy(N, x,1, y,1);
	
    zCreate_Dense_Matrix(&B, N, 1, (doublecomplex*)y, N, SLU_DN, SLU_Z, SLU_GE);
	zgstrs(NOTRANS, &superlu_data.L, &superlu_data.U,
		superlu_data.perm_c, superlu_data.perm_r,
		&B, &superlu_data.stat, &info);
	//std::cout << "info  = " << info << std::endl;
}

int SPB::BandSolver_Ez::MakeANumeric(){
	return 0;
}
int SPB::BandSolver_Ez::UpdateA(const double k[2]){
	//temp hacky
	InvalidateByStructure();
	if(NULL == ind){
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
