#include "SPB.hpp"
#include <cstdlib>
#include <limits>
#include <ctime>
#define RNP_OUTPUT_MATHEMATICA
//#define RNP_OUTPUT_MATLAB
#include <IO.h>
#define RNP_SPARSE_USE_IO
#include <Sparse.h>
#include <fftw3.h>
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
typedef RNP::Sparse::TCCSMatrix<complex_t> sparse_t;

#define DIVH_OFF (0*Ngrid)
#define EX_OFF   (undefined)
#define EY_OFF   (undefined)
#define EZ_OFF   (1*Ngrid)
#define HX_OFF   (2*Ngrid)
#define HY_OFF   (3*Ngrid)
#define HZ_OFF   (undefined)

#define IDX(i,j) (res[1]*(i)+(j))

class SPB::BandSolver_Ez::Impl{
	size_t N;
	int *ind;
	sparse_t *A, *B;
	friend class SPB::BandSolver_Ez;
	complex_t *eps_z_fft;
	bool structure_changed_since_last_solve;
public:
	Impl():N(0),ind(NULL),A(NULL),B(NULL),eps_z_fft(NULL),
		structure_changed_since_last_solve(true)
	{
	}
	~Impl(){
		free(ind);
		fftw_free(eps_z_fft);
		if(NULL != A){ delete A; }
		if(NULL != B){ delete B; }
	}
};

void SPB::BandSolver_Ez::ShiftInv(const complex_t &shift, const complex_t *x, complex_t *y) const{
	const int Ngrid = res[0]*res[1];
	size_t n = 4*Ngrid;
	
	RNP::TBLAS::Copy(n, x,1, y,1);
	
	// Invert V and P fields
	// For zero shift, the main E/H field block matrix inverse is unchanged
	// v = inv(D) h
	// u = inv(C - W^H inv(D) W) (g - W^H v)
	for(int i = 0; i < res[0]; ++i){
		for(int j = 0; j < res[1]; ++j){
			const int tag = impl->ind[2*IDX(i,j)+1];
			if(tag < 0){ continue; }
			const int row0 = impl->ind[2*IDX(i,j)+0];
			const Material &mat = material[tag];
			const int np = mat.poles.size();
			for(int p = 0; p < np; ++p){
				const LorentzPole &pole = mat.poles[p];
				const complex_t iwp = complex_t(0,pole.omega_p) * mat.eps_inf.value[8];
				const complex_t i_w0 = complex_t(0,1./pole.omega_0) / mat.eps_inf.value[8];
				const complex_t iG_w0w0 = i_w0 * pole.Gamma / pole.omega_0;
				y[row0 + 2*p + 0] = -i_w0 * x[row0 + 2*p + 1];
				y[row0 + 2*p + 1] = i_w0 * x[row0 + 2*p + 0] + iG_w0w0 * x[row0 + 2*p + 1];
				y[EZ_OFF + IDX(i,j)] += iwp * y[row0 + 2*p + 0];
			}
		}
	}
	
	// Data layout: divH, Ez, Hx, Hy
	fftw_plan plan_forward = fftw_plan_many_dft(
		2/*rank*/, res, 4 /*howmany*/,
		(fftw_complex*)y, NULL/*inembed*/,
		1/*istride*/, Ngrid/*idist*/,
		(fftw_complex*)y, NULL/*onembed*/,
		1/*ostride*/, Ngrid/*odist*/,
		FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_backward = fftw_plan_many_dft(
		2/*rank*/, res, 4 /*howmany*/,
		(fftw_complex*)y, NULL/*inembed*/,
		1/*istride*/, Ngrid/*idist*/,
		(fftw_complex*)y, NULL/*onembed*/,
		1/*ostride*/, Ngrid/*odist*/,
		FFTW_FORWARD, FFTW_ESTIMATE);
	const double kshiftsign = 1.0;
	for(int i = 0; i < res[0]; ++i){
		for(int j = 0; j < res[1]; ++j){
			double phase = kshiftsign*2*M_PI*(last_k[0]*(double)i/res[0] + last_k[1]*(double)j/res[1]);
			y[HX_OFF + IDX(i,j)] *= complex_t(cos(phase), sin(phase));
			y[HY_OFF + IDX(i,j)] *= complex_t(cos(phase), sin(phase));
			y[EZ_OFF + IDX(i,j)] *= complex_t(cos(phase), sin(phase));
			y[DIVH_OFF+IDX(i,j)] *= complex_t(cos(phase), sin(phase));
		}
	}
	fftw_execute(plan_forward);

	complex_t A[4*4], b[4];
	size_t ipiv[4];
	for(int i = 0; i < res[0]; ++i){
		const int fi = (i > res[0]/2 ? i-res[0] : i);
		for(int j = 0; j < res[1]; ++j){
			const int fj = (j > res[1]/2 ? j-res[1] : j);
			double kpG[2] = {
				(L.Lk[0]*(last_k[0]+fi) + L.Lk[2]*(last_k[1]+fj)),
				(L.Lk[1]*(last_k[0]+fi) + L.Lk[3]*(last_k[1]+fj))
			};
			kpG[0] *= 2*M_PI;
			kpG[1] *= 2*M_PI;
			const double klen2 = kpG[0]*kpG[0] + kpG[1]*kpG[1];
			const double klen = sqrt(klen2);
			
			// At the Gamma point, project out the constant basis vector
			if(klen < std::numeric_limits<double>::epsilon() * L.CharacteristicKLength()){
				y[DIVH_OFF+IDX(i,j)] = 0;
				y[EZ_OFF + IDX(i,j)] = 0;
				y[HX_OFF + IDX(i,j)] = 0;
				y[HY_OFF + IDX(i,j)] = 0;
				continue;
			}
			
			// [   0      0      k.        0         0     ] [dvH] = [dvH]
			// [   0   -q eps   -k x       0       -i wp   ] [ E ] = [ E ]
			// [   k     k x   -q mu       0         0     ] [ H ]   [ H ]
			// [   0      0      0      -q eta      i w0   ] [ P ]   [ P ]
			// [   0    i wp     0     -i w0 eta   -q eta  ] [ V ]   [ V ]
			//                                                       given
			memset(A, 0, sizeof(complex_t)*4*4);
			
// Forward and backward differences
#define FDIFF(VEC,D) ((std::exp(complex_t(0,-(VEC)[D]/res[D]))-1.) * (double)res[D])
#define BDIFF(VEC,D) ((1.-std::exp(complex_t(0,(VEC)[D]/res[D]))) * (double)res[D])
			static const complex_t I(0.,1.);
			
			A[2+1*4] = -I*FDIFF(kpG,1);
			A[2+0*4] = I*BDIFF(kpG,0);
			A[3+1*4] = I*FDIFF(kpG,0);
			A[3+0*4] = I*BDIFF(kpG,1);
			A[1+2*4] = -I*BDIFF(kpG,1);
			A[1+3*4] =  I*BDIFF(kpG,0);
			A[0+2*4] = I*FDIFF(kpG,0);
			A[0+3*4] = I*FDIFF(kpG,1);
			
			b[0] = y[DIVH_OFF+IDX(i,j)];
			b[1] = y[EZ_OFF + IDX(i,j)];
			b[2] = y[HX_OFF + IDX(i,j)];
			b[3] = y[HY_OFF + IDX(i,j)];
			
			RNP::LinearSolve<'N'>(4,1, A,4, b,4);
			
			y[DIVH_OFF+IDX(i,j)] = b[0] / ((double)Ngrid);
			y[EZ_OFF + IDX(i,j)] = b[1] / ((double)Ngrid);
			y[HX_OFF + IDX(i,j)] = b[2] / ((double)Ngrid);
			y[HY_OFF + IDX(i,j)] = b[3] / ((double)Ngrid);
		}
	}
	fftw_execute(plan_backward);
	fftw_destroy_plan(plan_forward);
	fftw_destroy_plan(plan_backward);
	
	for(int i = 0; i < res[0]; ++i){
		for(int j = 0; j < res[1]; ++j){
			double phase = -kshiftsign*2*M_PI*(last_k[0]*(double)i/res[0] + last_k[1]*(double)j/res[1]);
			y[DIVH_OFF+IDX(i,j)] *= complex_t(cos(phase), sin(phase));
			y[EZ_OFF + IDX(i,j)] *= complex_t(cos(phase), sin(phase));
			y[HX_OFF + IDX(i,j)] *= complex_t(cos(phase), sin(phase));
			y[HY_OFF + IDX(i,j)] *= complex_t(cos(phase), sin(phase));
		}
	}
	
	// v += inv(D) W u
	for(int i = 0; i < res[0]; ++i){
		for(int j = 0; j < res[1]; ++j){
			const int tag = impl->ind[2*IDX(i,j)+1];
			if(tag < 0){ continue; }
			const int row0 = impl->ind[2*IDX(i,j)+0];
			const Material &mat = material[tag];
			const int np = mat.poles.size();
			for(int p = 0; p < np; ++p){
				const LorentzPole &pole = mat.poles[p];
				y[row0 + 2*p + 0] -= (pole.omega_p/pole.omega_0) * y[EZ_OFF + IDX(i,j)];
			}
		}
	}
}

SPB::BandSolver_Ez::BandSolver_Ez(double Lr[4]):BandSolver(Lattice2(Lr)),L(Lr){
	impl = new Impl();
}

SPB::BandSolver_Ez::~BandSolver_Ez(){
	delete impl;
}

size_t SPB::BandSolver_Ez::GetSize() const{
	if(NULL != impl){
		return impl->N;
	}else{
		return 0;
	}
}
void SPB::BandSolver_Ez::Aop(const std::complex<double> *x, std::complex<double> *y) const{
	RNP::Sparse::MultMV<'N'>(*(impl->A), x, y);
}
void SPB::BandSolver_Ez::Bop(const std::complex<double> *x, std::complex<double> *y) const{
	RNP::Sparse::MultMV<'N'>(*(impl->B), x, y);
}
void SPB::BandSolver_Ez::PrintField(const std::complex<double> *y, const char *filename) const{
	const int Ngrid = res[0]*res[1];
	std::ostream *out = &std::cout;
	if(NULL != filename){
		out = new std::ofstream(filename);
	}
	
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
	}
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

	// Prepare the indexing
	const size_t Ngrid = res[0] * res[1];
	
	if(impl->structure_changed_since_last_solve){
		free(impl->ind);
		impl->ind = (int*)malloc(sizeof(int) * 2*Ngrid);

		fftw_free(impl->eps_z_fft);
		impl->eps_z_fft = (complex_t*)fftw_malloc(sizeof(complex_t)*Ngrid);
		
		size_t next_index = 0;
		for(int i = 0; i < res[0]; ++i){
			const double fi = ((double)i/(double)res[0]) - 0.5;
			for(int j = 0; j < res[1]; ++j){
				const double fj = ((double)j/(double)res[1]) - 0.5;
				impl->ind[2*IDX(i,j)+0] = 4*Ngrid+next_index;
				
				// get material of this cell (simple pointwise check)
				int tag, num_poles;
				//if(2 == dim){
					double p[2] = {
						L.Lr[0]*fi + L.Lr[2]*fj,
						L.Lr[1]*fi + L.Lr[3]*fj
					};
					if(!shapeset.QueryPt(p, &tag)){
						tag = -1;
					}
				/*}else{
					double p[3] = {
						L.Lr[0]*fi + L.Lr[3]*fj + L.Lr[6]*fk,
						L.Lr[1]*fi + L.Lr[4]*fj + L.Lr[7]*fk,
						L.Lr[2]*fi + L.Lr[5]*fj + L.Lr[8]*fk
					};
					if(ShapeSet3_query_pt(shapeset.d3, p, NULL, &tag)){
					}else{
						tag = -1;
					}
				}*/
				if(-1 == tag){
					num_poles = 0;
					impl->eps_z_fft[IDX(i,j)] = 1.;
				}else{
					num_poles = material[tag].poles.size();
					impl->eps_z_fft[IDX(i,j)] = material[tag].eps_inf.value[8];
std::cout << i << "\t" << j << "\t" << impl->eps_z_fft[IDX(i,j)] << "\t" << num_poles << std::endl;
				}
				impl->ind[2*IDX(i,j)+1] = tag;
				// update next index
				next_index += 2*num_poles;
			}
		}
		//impl->N = 4*Ngrid + 3*zero_constraint + next_index;
		impl->N = 4*Ngrid + next_index;


		/*
		switch(pol){
		case 1:
			// Hx,Hy,Ez, divH
			N = (3+1)*Ngrid + 3*zero_constraint + next_index;
			break;
		case 2:
			// Hz,Ex,Ey (Hz is already div-free)
			N = (3+0)*Ngrid + 3*zero_constraint + next_index;
			break;
		default:
			// Hx,Hy,Hz,Ex,Ey,Ez, divH
			N = (6+1)*Ngrid + 6*zero_constraint + next_index;
			break;
		}*/
		
		fftw_plan plan_eps = fftw_plan_many_dft(
			2/*rank*/, res, 1 /*howmany*/,
			(fftw_complex*)impl->eps_z_fft, NULL/*inembed*/,
			1/*istride*/, Ngrid/*idist*/,
			(fftw_complex*)impl->eps_z_fft, NULL/*onembed*/,
			1/*ostride*/, Ngrid/*odist*/,
			FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan_eps);
		fftw_destroy_plan(plan_eps);
		impl->structure_changed_since_last_solve = false;
	}
	
	sparse_t::entry_map_t Amap;
	sparse_t::entry_map_t Bmap;
	
	{
		const double Lrl[2] = {
			hypot(L.Lr[0], L.Lr[1]),
			hypot(L.Lr[2], L.Lr[3])
		};
		const double idr[2] = {
			(double)res[0] / Lrl[0],
			(double)res[1] / Lrl[1]
		};
				
		const complex_t Bloch[2] = {
			complex_t(cos(k[0]*2*M_PI), sin(k[0]*2*M_PI)),
			complex_t(cos(k[1]*2*M_PI), sin(k[1]*2*M_PI))
		};
		for(int i = 0; i < res[0]; ++i){
			for(int j = 0; j < res[1]; ++j){
				size_t row, col;
				complex_t coeff;
				
				const int curmat = impl->ind[2*IDX(i,j)+1];
				complex_t eps_z(1.);
				if(curmat >= 0){
					eps_z = material[curmat].eps_inf.value[8];
				}
				
#define ASET(ROW,COL,COEFF) Amap[sparse_t::index_t((ROW),(COL))] = (COEFF)
#define BSET(ROW,COL,COEFF) Bmap[sparse_t::index_t((ROW),(COL))] = (COEFF)
				// divH ~ dx Hx + dy Hy + dz Hz
				// E ~ -i wp V
				// V ~ +i wp E - i G V - i w0 P
				// P ~ +i w0 V
			
				//for(size_t idbg=0;idbg<ne+nh+1;++idbg){
					//ASET(row0+idbg,row0+idbg,1); // for debugging
				//}
			
				// Hx ~ -i dy Ez
				// Hy ~ +i dx Ez
				// Ez ~ -i dy Hx + i dx Hy

				// Hx = complex_t(0,-idr[1]) * (Ez[i,j+1,k] - Ez[i,j,k])
				row = HX_OFF + IDX(i,j);
				coeff = complex_t(0,-idr[1]);
				col = EZ_OFF + IDX(i,j); // Ez
				ASET(row,col, -coeff);
				if(j+1 == res[1]){
					col = EZ_OFF + IDX(i,0); // Ez
					ASET(row,col, coeff/Bloch[1]);
				}else{
					col = EZ_OFF + IDX(i,j+1); // Ez
					ASET(row,col, coeff);
				}
				BSET(row,row, 1);
				
				// Hy = complex_t(0, idr[0]) * (Ez[i+1,j,k] - Ez[i,j,k])
				row = HY_OFF + IDX(i,j);
				coeff = complex_t(0, idr[0]);
				col = EZ_OFF + IDX(i,j); // Ez
				ASET(row,col, -coeff);
				if(i+1 == res[0]){
					col = EZ_OFF + IDX(0,j); // Ez
					ASET(row,col, coeff/Bloch[0]);
				}else{
					col = EZ_OFF + IDX(i+1,j); // Ez
					ASET(row,col, coeff);
				}
				BSET(row,row, 1);
				
				// divH = idr[0] * (Hx[i+1,j,k] - Hx[i,j,k])
				//      + idr[1] * (Hy[i,j+1,k] - Hx[i,j,k])
				row = DIVH_OFF + IDX(i,j);
				coeff = complex_t(0,idr[0]);
				col = HX_OFF + IDX(i,j); // Hx
				ASET(row,col, -coeff);
				ASET(col,row, -std::conj(coeff));
				if(i+1 == res[0]){
					col = HX_OFF + IDX(0,j); // Hx
					ASET(row,col, coeff/Bloch[0]);
					ASET(col,row, std::conj(coeff/Bloch[0]));
				}else{
					col = HX_OFF + IDX(i+1,j); // Hx
					ASET(row,col, coeff);
					ASET(col,row, std::conj(coeff));
				}
				
				coeff = complex_t(0,idr[1]);
				col = HY_OFF + IDX(i,j); // Hy
				ASET(row,col, -coeff);
				ASET(col,row, -std::conj(coeff));
				if(j+1 == res[1]){
					col = HY_OFF + IDX(i,0); // Hy
					ASET(row,col, coeff/Bloch[1]);
					ASET(col,row, std::conj(coeff/Bloch[1]));
				}else{
					col = HY_OFF + IDX(i,j+1); // Hy
					ASET(row,col, coeff);
					ASET(col,row, std::conj(coeff));
				}
				BSET(row,row, 0);

				// Ez = complex_t(0,-idr[1]) * (Hx[i,j,k] - Hx[i,j-1,k])
				//    + complex_t(0, idr[0]) * (Hy[i,j,k] - Hy[i-1,j,k])
				row = EZ_OFF + IDX(i,j);
				
				coeff = complex_t(0,-idr[1]);
				col = HX_OFF + IDX(i,j); // Hx
				ASET(row,col, coeff);
				if(0 == j){
					col = HX_OFF + IDX(i,res[1]-1); // Hx
					ASET(row,col, -coeff*Bloch[1]);
				}else{
					col = HX_OFF + IDX(i,j-1); // Hx
					ASET(row,col, -coeff);
				}
				
				coeff = complex_t(0, idr[0]);
				col = HY_OFF + IDX(i,j); // Hy
				ASET(row,col, coeff);
				if(0 == i){
					col = HY_OFF + IDX(res[0]-1,j); // Hy
					ASET(row,col, -coeff*Bloch[0]);
				}else{
					col = HY_OFF + IDX(i-1,j); // Hy
					ASET(row,col, -coeff);
				}
				BSET(row,row, eps_z);
				
				if(curmat >= 0){
					const int row0 = impl->ind[2*IDX(i,j)+0];
					const Material &m = material[curmat];
					const size_t np = m.poles.size();
					for(size_t p = 0; p < np; ++p){
						row = row0 + 2*p + 0; // V_p
						coeff = complex_t(0, m.poles[p].omega_p) * eps_z;
						col = EZ_OFF + IDX(i,j); // E
						ASET(row,col, coeff);
						ASET(col,row, std::conj(coeff));
						
						if(0 != m.poles[p].Gamma){
							coeff = complex_t(0,-m.poles[p].Gamma) * eps_z;
							ASET(row,row, coeff);
						}
						BSET(row,row, 1);
						
						
						coeff = complex_t(0, -m.poles[p].omega_0) * eps_z;
						col = row0 + 2*p + 1; // P
						ASET(row,col, coeff);
						ASET(col,row, std::conj(coeff));
						BSET(col,col, 1);
					}
				}
				
				/*
				}else if(2 == pol){
					// Hz ~ +i dy Ex - i dx Ey
					// Ex ~ +i dy Hz
					// Ey ~ -i dx Hz
				}else{
					// Hx ~ +i dz Ey - i dy Ez
					// Hy ~ -i dz Ex + i dx Ez
					// Hz ~ +i dy Ex - i dx Ey
					// Ex ~ -i dz Hy + i dy Hz
					// Ey ~ +i dz Hx - i dx Hz
					// Ez ~ -i dy Hx + i dx Hy
				}*/
			}
		}
	}
	impl->A = new sparse_t(impl->N,impl->N, Amap);
	impl->B = new sparse_t(impl->N,impl->N, Bmap);
	
	if(0){
		std::cout << "A="; RNP::Sparse::PrintSparseMatrix(*(impl->A)) << ";" << std::endl;
		std::cout << "B="; RNP::Sparse::PrintSparseMatrix(*(impl->B)) << ";" << std::endl;
		exit(0);
	}
	/*
	complex_t *tmp = new complex_t[4*Ngrid];
	complex_t *tmp2 = new complex_t[16*Ngrid*Ngrid];
	for(size_t i = 0; i < res[0]; ++i){
		for(size_t j = 0; j < res[1]; ++j){
			tmp[IDX(i,j)] = 0;
		}
	}
	for(size_t i = 0; i < res[0]; ++i){
		for(size_t j = 0; j < res[1]; ++j){
			tmp[IDX(i,j)] = 1;
			Precond(tmp, &tmp2[0+IDX(i,j)*Ngrid]);
			tmp[IDX(i,j)] = 0;
		}
	}
	delete [] tmp2;
	delete [] tmp;
	*/
	return solver->Solve();
	/*
	{
		const size_t n = 4*Ngrid;
		complex_t *x = (complex_t*)fftw_malloc(sizeof(complex_t)*n);
		complex_t *y = (complex_t*)fftw_malloc(sizeof(complex_t)*n);
		complex_t *z = (complex_t*)fftw_malloc(sizeof(complex_t)*n);
		const double theta = 0.6;
		
		memset(x, 0, sizeof(complex_t)*n);
		for(int i = 0; i < n; ++i){
			x[i] = frand();
		}
		std::cout << "x = "; RNP::IO::PrintVector(n, x, 1) << std::endl;
		
		Aop(x, y);
		Bop(x, z);
		RNP::TBLAS::Axpy(n, -theta, z,1, y,1);
		// At this point y = A*x-theta*B*x
		std::cout << "y = "; RNP::IO::PrintVector(n, y, 1) << std::endl;
		
		Op(n, theta, y, z);
		// At this point z should be the same as x
		
		std::cout << "z = "; RNP::IO::PrintVector(n, z, 1) << std::endl;
		
		RNP::TBLAS::Axpy(n, -1., x,1,z,1);
		std::cout << "diff = "; RNP::IO::PrintVector(n, z, 1) << std::endl;
		
		fftw_free(z);
		fftw_free(y);
		fftw_free(x);
	}*/
	
	
	/*
	size_t n_wanted = 10;
	size_t ncv = 2*n_wanted+1;
	SPB::complex_t *w = new SPB::complex_t[n_wanted+ncv*4*Ngrid];
	SPB::complex_t *v = w+n_wanted;
	int nconv = RNP::IRA::ShiftInvert(
		4*Ngrid, 0.0, &op_, &bv_,
		n_wanted, ncv, &RNP::IRA::LargestMagnitude,
		w, v, 4*Ngrid,
		NULL,
		NULL,
		(void*)this,
		(void*)this);
	for(size_t i = 0; i < n_wanted;++i){
		std::cout << w[i] << std::endl;
	}
	*/
}

// This is an operator that applies (A - shift*B) on a vector using the FFT
void SPB::BandSolver_Ez::OpForw(size_t n, const complex_t &shift, const complex_t *x, complex_t *y) const{
	const int Ngrid = res[0]*res[1];
	complex_t *t = (complex_t*)fftw_malloc(sizeof(complex_t)*n);
	// Data layout: Hx, Hy, Ez, divH
	fftw_plan plan_forward = fftw_plan_many_dft(
		2/*rank*/, res, 4 /*howmany*/,
		(fftw_complex*)t, NULL/*inembed*/,
		1/*istride*/, Ngrid/*idist*/,
		(fftw_complex*)t, NULL/*onembed*/,
		1/*ostride*/, Ngrid/*odist*/,
		FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_backward = fftw_plan_many_dft(
		2/*rank*/, res, 4 /*howmany*/,
		(fftw_complex*)t, NULL/*inembed*/,
		1/*istride*/, Ngrid/*idist*/,
		(fftw_complex*)t, NULL/*onembed*/,
		1/*ostride*/, Ngrid/*odist*/,
		FFTW_FORWARD, FFTW_ESTIMATE);
	
	const double kshiftsign = 1.0;
	
	RNP::TBLAS::Copy(n, x,1, t,1);
	for(int i = 0; i < res[0]; ++i){
		for(int j = 0; j < res[1]; ++j){
			double phase = kshiftsign*2*M_PI*(last_k[0]*((double)i)/res[0] + last_k[1]*((double)j)/res[1]);
			t[HX_OFF + IDX(i,j)] *= complex_t(cos(phase), sin(phase));
			t[HY_OFF + IDX(i,j)] *= complex_t(cos(phase), sin(phase));
			t[EZ_OFF + IDX(i,j)] *= complex_t(cos(phase), sin(phase));
			t[DIVH_OFF+IDX(i,j)] *= complex_t(cos(phase), sin(phase));
		}
	}
	fftw_execute(plan_forward);

	for(int i = 0; i < res[0]; ++i){
		const int fi = (i > res[0]/2 ? i-res[0] : i);
		for(int j = 0; j < res[1]; ++j){
			const int fj = (j > res[1]/2 ? j-res[1] : j);
			double kpG[2] = {
				(L.Lk[0]*(last_k[0]+fi) + L.Lk[2]*(last_k[1]+fj)),
				(L.Lk[1]*(last_k[0]+fi) + L.Lk[3]*(last_k[1]+fj))
			};
			kpG[0] *= 2*M_PI;
			kpG[1] *= 2*M_PI;
			const double klen2 = kpG[0]*kpG[0] + kpG[1]*kpG[1];
			const double klen = sqrt(klen2);
			if(klen < std::numeric_limits<double>::epsilon() * L.CharacteristicKLength()){ continue; }
			
			// [ -q mu     k x     k       0         0     ] [ H ]   [ H ]
			// [  -k x   -q eps    0       0       -i wp   ] [ E ] = [ E ]
			// [   k.       0      0       0         0     ] [dvH] = [dvH]
			// [   0        0      0    -q eta      i w0   ] [ P ]   [ P ]
			// [   0      i wp     0   -i w0 eta   -q eta  ] [ V ]   [ V ]
			//                                               given
/*
// Forward and backward differences
#define FDIFF(VEC,D) ((std::exp(complex_t(0,-(VEC)[D]/res[D]))-1.) * (double)res[D])
#define BDIFF(VEC,D) ((1.-std::exp(complex_t(0,(VEC)[D]/res[D]))) * (double)res[D])
			static const complex_t I(0.,1.);
			y[HX_OFF + IDX(i,j)] = -I*FDIFF(kpG,1)*t[EZ_OFF + IDX(i,j)] + I*BDIFF(kpG,0)*t[DIVH_OFF+IDX(i,j)];
			y[HX_OFF + IDX(i,j)] /= ((double)Ngrid);
			y[HY_OFF + IDX(i,j)] =  I*FDIFF(kpG,0)*t[EZ_OFF + IDX(i,j)] + I*BDIFF(kpG,1)*t[DIVH_OFF+IDX(i,j)];
			y[HY_OFF + IDX(i,j)] /= ((double)Ngrid);
			y[EZ_OFF + IDX(i,j)] = -I*BDIFF(kpG,1)*t[HX_OFF + IDX(i,j)] + I*BDIFF(kpG,0)*t[HY_OFF + IDX(i,j)];
			y[EZ_OFF + IDX(i,j)] /= ((double)Ngrid);
			y[DIVH_OFF+IDX(i,j)] =  I*FDIFF(kpG,0)*t[HX_OFF + IDX(i,j)] + I*FDIFF(kpG,1)*t[HY_OFF + IDX(i,j)];
			y[DIVH_OFF+IDX(i,j)] /= ((double)Ngrid);
			
			y[HX_OFF + IDX(i,j)] -= shift*t[HX_OFF + IDX(i,j)]/((double)Ngrid);
			y[HY_OFF + IDX(i,j)] -= shift*t[HY_OFF + IDX(i,j)]/((double)Ngrid);
			y[EZ_OFF + IDX(i,j)] -= shift*t[EZ_OFF + IDX(i,j)]/((double)Ngrid);
*/
			static const complex_t I(0.,1.);
			const size_t n_res = 0;
			const size_t nh = 4+2*n_res;
			complex_t *A = new complex_t[nh*nh+2*nh];
			complex_t *b = A+nh*nh;
			complex_t *c = b+nh;
			memset(A, 0, sizeof(complex_t)*nh*nh);
			
			A[0+2*nh] = -I*FDIFF(kpG,1);
			A[0+3*nh] = I*BDIFF(kpG,0);
			A[1+2*nh] = I*FDIFF(kpG,0);
			A[1+3*nh] = I*BDIFF(kpG,1);
			A[2+0*nh] = -I*BDIFF(kpG,1);
			A[2+1*nh] =  I*BDIFF(kpG,0);
			A[3+0*nh] = I*FDIFF(kpG,0);
			A[3+1*nh] = I*FDIFF(kpG,1);
			
			A[0+0*nh] = -shift;
			A[1+1*nh] = -shift;
			A[2+2*nh] = -shift;
			
			b[0] = t[HX_OFF + IDX(i,j)];
			b[1] = t[HY_OFF + IDX(i,j)];
			b[2] = t[EZ_OFF + IDX(i,j)];
			b[3] = t[DIVH_OFF+IDX(i,j)];
			
			RNP::TBLAS::MultMV<'N'>(nh,nh, 1.,A,nh, b,1, 0.,c,1);
			
			y[HX_OFF + IDX(i,j)] = c[0] / ((double)Ngrid);
			y[HY_OFF + IDX(i,j)] = c[1] / ((double)Ngrid);
			y[EZ_OFF + IDX(i,j)] = c[2] / ((double)Ngrid);
			y[DIVH_OFF+IDX(i,j)] = c[3] / ((double)Ngrid);
			
			delete [] A;
		}
	}
	RNP::TBLAS::Copy(n, y,1, t,1);
	fftw_execute(plan_backward);
	fftw_destroy_plan(plan_forward);
	fftw_destroy_plan(plan_backward);
	RNP::TBLAS::Copy(n, t,1, y,1);
	
	for(int i = 0; i < res[0]; ++i){
		for(int j = 0; j < res[1]; ++j){
			double phase = -kshiftsign*2*M_PI*(last_k[0]*((double)i)/res[0] + last_k[1]*((double)j)/res[1]);
			y[HX_OFF + IDX(i,j)] *= complex_t(cos(phase), sin(phase));
			y[HY_OFF + IDX(i,j)] *= complex_t(cos(phase), sin(phase));
			y[EZ_OFF + IDX(i,j)] *= complex_t(cos(phase), sin(phase));
			y[DIVH_OFF+IDX(i,j)] *= complex_t(cos(phase), sin(phase));
			/*
			y[HX_OFF + IDX(i,j)] -= shift*x[HX_OFF + IDX(i,j)];
			y[HY_OFF + IDX(i,j)] -= shift*x[HY_OFF + IDX(i,j)];
			y[EZ_OFF + IDX(i,j)] -= shift*x[EZ_OFF + IDX(i,j)];*/
		}
	}
	fftw_free(t);
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
