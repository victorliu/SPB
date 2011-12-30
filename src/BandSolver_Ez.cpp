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
				//ASET(col,col, -assembly_data.shift * Bvals);
				
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
	N(0),ldl(),
	valid_A_numeric(false)
{
	cell2ind = NULL;
	ind2cell = NULL;
	matind = NULL;
	npoles = NULL;
	assembly_data.shift = 2*M_PI*0.2;
}

SPB::BandSolver_Ez::~BandSolver_Ez(){
	free(cell2ind);
	free(ind2cell);
	free(matind);
	free(npoles);
}

size_t SPB::BandSolver_Ez::GetSize() const{
	return N;
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
	
	return 0;
}


int SPB::BandSolver_Ez::GetMaxBlockSize() const{
	return assembly_data.max_block_size;
}
int SPB::BandSolver_Ez::GetMaxNNZPerRow() const{
	return assembly_data.max_nnz_per_row;
}
int SPB::BandSolver_Ez::BeginBlockSymbolic() const{
	assembly_data.p = 0;
	return res[0]*res[1];
}
int SPB::BandSolver_Ez::BeginBlockNumeric() const{
	assembly_data.p = 0;
	return res[0]*res[1];
}
int SPB::BandSolver_Ez::GetNextBlockSymbolic(int *rowptr, int *colind) const{
	return GetNextBlockNumeric(rowptr, colind, NULL);
}
int SPB::BandSolver_Ez::GetNextBlockNumeric(int *rowptr, int *colind, complex_t *value) const{
	size_t row;
	complex_t coeff;
	int next_col = 0, row_count = 0;

	if(assembly_data.p >= res[0]*res[1]){ return 0; }
	int q = ind2cell[2*assembly_data.p+1];
	int i,j;
	UNIDX(q,i,j);
	
//std::cerr << "  Getting block for p,q,i,j=" << assembly_data.p << "," << q << "," << i << "," << j << std::endl;

#define NEWROW() do{ *rowptr = next_col; ++rowptr; ++row_count; }while(0)
#define SETCOL(I,J,W) do{ \
		*colind = cell2ind[IDX((I),(J))]+(W); colind++; next_col++; \
	}while(0)
#define SETVAL(VAL) do{ \
		if(NULL != value){ *value = (VAL); value++; } \
	}while(0)

	const int curmat = matind[q]-2;
	double eps_z = epsval[epsind[q]];
	const double Lrl[2] = {
		hypot(L.Lr[0], L.Lr[1]),
		hypot(L.Lr[2], L.Lr[3])
	};
	const double idr[2] = {
		(double)res[0] / Lrl[0],
		(double)res[1] / Lrl[1]
	};

	
	// Hy = idr[0] * (Ez[i+1,j,k] - Ez[i,j,k])
	//    + idr[1] * (divH[i,j,k] - divH[i,j-1,k])
	NEWROW();
	{
		// diagonal
		SETCOL(i, j, HY_OFF);
		SETVAL(-assembly_data.shift);
		// Ez coupling
		coeff = complex_t(0,idr[0]);
		if(i+1 == res[0]){
			SETCOL(0, j, EZ_OFF);
			SETVAL(coeff*assembly_data.Bloch[0]);
		}else{
			SETCOL(i+1, j, EZ_OFF);
			SETVAL(coeff);
		}
		SETCOL(i, j, EZ_OFF);
		SETVAL(-coeff);
		// divH coupling
		coeff = complex_t(0,idr[1]);
		if(0 == j){
			SETCOL(i, res[1]-1, DIVH_OFF);
			SETVAL(-coeff/assembly_data.Bloch[1]);
		}else{
			SETCOL(i, j-1, DIVH_OFF);
			SETVAL(-coeff);
		}
		SETCOL(i, j, DIVH_OFF);
		SETVAL(coeff);
	}
	// Ez = -idr[1] * (Hx[i,j,k] - Hx[i,j-1,k])
	//    +  idr[0] * (Hy[i,j,k] - Hy[i-1,j,k])
	NEWROW();
	{
		// diagonal
		SETCOL(i, j, EZ_OFF);
		SETVAL(-assembly_data.shift*eps_z);
		// Hy coupling
		coeff = complex_t(0,idr[0]);
		if(0 == i){
			SETCOL(res[0]-1, j, HY_OFF);
			SETVAL(-coeff/assembly_data.Bloch[0]);
		}else{
			SETCOL(i-1, j, HY_OFF);
			SETVAL(-coeff);
		}
		SETCOL(i, j, HY_OFF);
		SETVAL(coeff);
		// Hx coupling
		coeff = complex_t(0,idr[1]);
		if(0 == j){
			SETCOL(i, res[1]-1, HX_OFF);
			SETVAL(coeff/assembly_data.Bloch[1]);
		}else{
			SETCOL(i, j-1, HX_OFF);
			SETVAL(coeff);
		}
		SETCOL(i, j, HX_OFF);
		SETVAL(-coeff);
	}
	// Hx += -idr[1] * (Ez[i,j+1,k] - Ez[i,j,k])
	//    +   idr[0] * (divH[i,j,k] - divH[i-1,j,k])
	NEWROW();
	{
		// diagonal
		SETCOL(i, j, HX_OFF);
		SETVAL(-assembly_data.shift);
		// Ez coupling
		coeff = complex_t(0,idr[1]);
		if(j+1 == res[1]){
			SETCOL(i, 0, EZ_OFF);
			SETVAL(-coeff*assembly_data.Bloch[1]);
		}else{
			SETCOL(i, j+1, EZ_OFF);
			SETVAL(-coeff);
		}
		SETCOL(i, j, EZ_OFF);
		SETVAL(coeff);
		// divH coupling
		coeff = complex_t(0,idr[0]);
		if(0 == i){
			SETCOL(res[0]-1, j, DIVH_OFF);
			SETVAL(-coeff/assembly_data.Bloch[0]);
		}else{
			SETCOL(i-1, j, DIVH_OFF);
			SETVAL(-coeff);
		}
		SETCOL(i, j, DIVH_OFF);
		SETVAL(coeff);
	}
	// divH = idr[0] * (Hx[i+1,j,k] - Hx[i,j,k])
	//      + idr[1] * (Hy[i,j+1,k] - Hy[i,j,k])
	NEWROW();
	{
		// diagonal
		SETCOL(i, j, DIVH_OFF);
		SETVAL(0.);
		// Hy coupling
		coeff = complex_t(0,idr[1]);
		if(j+1 == res[0]){
			SETCOL(i, 0, HY_OFF);
			SETVAL(coeff*assembly_data.Bloch[1]);
		}else{
			SETCOL(i, j+1, HY_OFF);
			SETVAL(coeff);
		}
		SETCOL(i, j, HY_OFF);
		SETVAL(-coeff);
		// Hx coupling
		coeff = complex_t(0,idr[0]);
		if(i+1 == res[0]){
			SETCOL(0, j, HX_OFF);
			SETVAL(coeff*assembly_data.Bloch[0]);
		}else{
			SETCOL(i+1, j, HX_OFF);
			SETVAL(coeff);
		}
		SETCOL(i, j, HX_OFF);
		SETVAL(-coeff);
	}
	
	int matbits = matind[q];
	while(matbits & 0xF){
		int m = (matbits & 0xF)-1;
		if(m > 0){
			m--;
			const Material& curmat = material[m];
			int col = NUM_EH;
			for(int pi = 0; pi < curmat.poles.size(); ++pi){
				const LorentzPole &pole = material[m].poles[pi];
				// P = -i w0 eps V
				NEWROW();
				{
					SETCOL(i, j, col);
					SETVAL(-assembly_data.shift*eps_z);
					SETCOL(i, j, col+1);
					SETVAL(complex_t(0.,-pole.omega_0) * eps_z);
				}
				// V = i w0 eps P - i wp eps E
				NEWROW();
				{
					SETCOL(i, j, col+1);
					SETVAL(-assembly_data.shift*eps_z);
					SETCOL(i, j, col);
					SETVAL(complex_t(0., pole.omega_0) * eps_z);
					SETCOL(i, j, EZ_OFF);
					SETVAL(complex_t(0.,-pole.omega_p) * eps_z);
				}
				col += 2;
			}
		}
		matbits >>= 4;
	}
	// finalize with NNZ
	*rowptr = next_col;
	
	assembly_data.p++;
	return row_count;
}

int SPB::BandSolver_Ez::MakeASymbolic(){
	const size_t Ngrid = res[0] * res[1];

	if(NULL != ind2cell){ free(ind2cell); }
	if(NULL != cell2ind){ free(cell2ind); }
	if(NULL != matind){ free(matind); }
	if(NULL != npoles){ free(npoles); }
	ind2cell = (int*)malloc(sizeof(int) * 2*Ngrid);
	cell2ind = (int*)malloc(sizeof(int) * Ngrid);
	matind = (int*)malloc(sizeof(int) * Ngrid);
	npoles = (int*)malloc(sizeof(int) * Ngrid);
	
	const double Lr[4] = {
		L.Lr[0], L.Lr[1],
		L.Lr[2], L.Lr[3]
	};

	// Build mesh
	//   u edge
	mesh.edge[4*0+2*0+0] = 0;
	mesh.edge[4*0+2*0+1] = 0;
	mesh.edge[4*0+2*1+0] = 1;
	mesh.edge[4*0+2*1+1] = 0;
	//   v edge
	mesh.edge[4*1+2*0+0] = 0;
	mesh.edge[4*1+2*0+1] = 0;
	mesh.edge[4*1+2*1+0] = 0;
	mesh.edge[4*1+2*1+1] = 1;
	const double uv = Lr[0]*Lr[2] + Lr[1]*Lr[3];
	const double uxv = Lr[0]*Lr[3] - Lr[1]*Lr[2];
	const double ulen = hypot(Lr[0], Lr[1]);
	const double vlen = hypot(Lr[2], Lr[3]);
	if(fabs(uv) < std::numeric_limits<double>::epsilon()*uxv){
		mesh.star_mu[0] = ulen/vlen;
		mesh.star_mu[1] = vlen/ulen;
	}else{
		mesh.n_edges = 3;
		double uvw[2], w[2];
		if(uv < 0){ // wide angle between u and v
			// edge w is u+v
			w[0] = Lr[0] + Lr[2];
			w[1] = Lr[1] + Lr[3];
			mesh.edge[4*2+2*0+0] = 0;
			mesh.edge[4*2+2*0+1] = 0;
			mesh.edge[4*2+2*1+0] = 1;
			mesh.edge[4*2+2*1+1] = 1;
			// star_mu[0] = ulen/dual_ulen
			// dual_ulen = circum(0,u,u+v) - circum(0,-v,u)
			// where circum(0,a,b) = ((|a|^2 b - |b|^2 a) x (axb))/(2 |axb|^2)
			// this simplifies to
			// dual_ulen = ||u|^2(v+w)-(w^2+v^2)u| / (2|uxv|)
			// generally, this is
			// dual_*len = (|*|^2(u+v+w) - (v^2+u^2+w^2)*) / (2 |uxv|)
			uvw[0] = 2*w[0];
			uvw[1] = 2*w[1];
		}else{
			// edge w is v-u
			mesh.edge[4*2+2*0+0] = 1;
			mesh.edge[4*2+2*0+1] = 0;
			mesh.edge[4*2+2*1+0] = 0;
			mesh.edge[4*2+2*1+1] = 1;
			uvw[0] = 2*Lr[2];
			uvw[1] = 2*Lr[3];
		}
		const double wlen = hypot(w[0],w[1]);
		const double uvw2 = ulen*ulen + vlen*vlen + wlen*wlen;
		mesh.star_mu[0] = 2*uxv*ulen / hypot(ulen*ulen*uvw[0] - uvw2*Lr[0], ulen*ulen*uvw[1] - uvw2*Lr[1]);
		mesh.star_mu[1] = 2*uxv*vlen / hypot(vlen*vlen*uvw[0] - uvw2*Lr[2], vlen*vlen*uvw[1] - uvw2*Lr[3]);
		mesh.star_mu[2] = 2*uxv*wlen / hypot(wlen*wlen*uvw[0] - uvw2* w[0], wlen*wlen*uvw[1] - uvw2* w[1]);
	}
	
	{
		double use_k[2] = { last_k[0], last_k[1] };
		assembly_data.Bloch[0] = complex_t(cos(use_k[0]*2*M_PI), sin(use_k[0]*2*M_PI));
		assembly_data.Bloch[1] = complex_t(cos(use_k[1]*2*M_PI), sin(use_k[1]*2*M_PI));
	}

	// Prepare the indexing
	epsval.clear();
	epsind.resize(Ngrid);
	{
		std::map<double,int> epsmap;
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
				//ind2cell[2*q+1] = q;
				
				// get materials of this cell (for now, simple pointwise check)
				matind[q] = 0;
				int tag, num_poles = 0;
				double p[2] = {
					Lr[0]*fi + Lr[2]*fj,
					Lr[1]*fi + Lr[3]*fj
				};
				if(!shapeset.QueryPt(p, &tag)){
					tag = -1;
				}
				
				double eps_z = 1.;
				if(tag >= 0){
					const Material& curmat = material[tag];
					for(int pi = 0; pi < curmat.poles.size(); ++pi){
						num_poles++;
						num_poles++;
					}
					eps_z = curmat.eps_inf.value[8].real();
				}
				npoles[q] = num_poles;
				if(tag <= 13){
					matind[q] = tag+2;
				}
				
				std::map<double,int>::const_iterator mapiter = epsmap.find(eps_z);
				if(mapiter != epsmap.end()){
					epsind[q] = mapiter->second;
				}else{
					epsind[q] = epsval.size();
					epsval.push_back(eps_z);
					epsmap[eps_z] = epsind[q];
				}
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
	
	assembly_data.max_nnz_per_row = 16;
	assembly_data.max_block_size = NUM_EH + 2*4; // arbitrary hack

	ldl.Factorize(*this);
	
	valid_A_numeric = true;
	
	return 0;
}

void SPB::BandSolver_Ez::ShiftInv(const complex_t &shift, const complex_t *x, complex_t *y) const{
	int info;
	RNP::TBLAS::Copy(N, x,1, y,1);
	ldl.Solve(y);
}
void SPB::BandSolver_Ez::Bop(const std::complex<double> *x, std::complex<double> *y) const{
	for(int i = 0; i < res[0]; ++i){
		for(int j = 0; j < res[1]; ++j){
			const int q = IDX(i,j);
			const int row0 = cell2ind[q];
			const int curmat = matind[q]-2;
			double eps_z = epsval[epsind[q]];
			y[row0+HY_OFF] = x[row0+HY_OFF];
			y[row0+EZ_OFF] = eps_z*x[row0+EZ_OFF];
			y[row0+HX_OFF] = x[row0+HX_OFF];
			y[row0+DIVH_OFF] = 0;
			for(int k = 0; k < npoles[q]; ++k){
				y[row0+NUM_EH+k] = eps_z*x[row0+NUM_EH+k];
			}
		}
	}
}

int SPB::BandSolver_Ez::UpdateA(const double k[2]){
	//temp hacky
	InvalidateByStructure();
	if(NULL == cell2ind){
		MakeASymbolic();
	}
	if(last_k[0] == k[0] && last_k[1] == k[1]){
		return 0;
	}
	// do update
	return 0;
}
void SPB::BandSolver_Ez::SetShift(double shift){
	assembly_data.shift = shift;
}

size_t SPB::BandSolver_Ez::GetProblemSize() const{
	return N;
}
int SPB::BandSolver_Ez::GetN() const{
	return N;
}
