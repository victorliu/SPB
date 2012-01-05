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
#include "libumesh.h"
}
#include "match.h"
extern "C"{
#include "mem.h"
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

#define HV_OFF   0
#define EZ_OFF   1
#define HU_OFF   2
#define DIVH_OFF 3
#define HW_OFF    4
#define DIVH2_OFF 5

#define IDX(i,j) (res[1]*(i)+(j))
#define UNIDX(q,i,j) do{ \
	(i) = (q) / res[1]; \
	(j) = (q) % res[1]; \
	}while(0)

SPB::BandSolver_Ez::BandSolver_Ez(double Lr[4]):BandSolver(Lattice2(Lr)),L(Lr),
	N(0)
{
	cell2ind = NULL;
	ind2cell = NULL;
	matind = NULL;
	npoles = NULL;
}

SPB::BandSolver_Ez::~BandSolver_Ez(){
	SPB_Free(npoles, "npoles at " __FILE__ ":" STR(__LINE__));
	SPB_Free(cell2ind, "cell2ind at " __FILE__ ":" STR(__LINE__));
	SPB_Free(ind2cell, "ind2cell at " __FILE__ ":" STR(__LINE__));
	SPB_Free(matind, "matind at " __FILE__ ":" STR(__LINE__));
}

void SPB::BandSolver_Ez::PrintField(const std::complex<double> *y, const char *filename) const{
	std::ostream *out = &std::cout;
	if(NULL != filename){
		out = new std::ofstream(filename);
	}
	for(int i = 0; i < res[0]; ++i){
		for(int j = 0; j < res[1]; ++j){
			int q = cell2ind[IDX(i,j)];
			(*out)
				<< i << "\t" << j
				 << "\t" << y[q+HU_OFF].real()
				 << "\t" << y[q+HU_OFF].imag()
				 << "\t" << y[q+HV_OFF].real()
				 << "\t" << y[q+HV_OFF].imag()
				 << "\t" << y[q+EZ_OFF].real()
				 << "\t" << y[q+EZ_OFF].imag()
				 << "\t" << y[q+DIVH_OFF].real()
				 << "\t" << y[q+DIVH_OFF].imag()
				 << "\n";
		}
		(*out) << "\n";
	}
	if(NULL != filename){
		delete out;
	}
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
	//SPB_VERB(9, ">GetNextBlockSymbolic\n");
	return GetNextBlockNumeric(rowptr, colind, NULL);
	//SPB_VERB(9, "<GetNextBlockSymbolic\n");
}
int SPB::BandSolver_Ez::GetNextBlockNumeric(int *rowptr, int *colind, complex_t *value) const{
	//SPB_VERB(9, ">GetNextBlockNumeric\n");
	const complex_t Im1(0.,1.);
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

	double eps_z = epsval[epsind[q]]*mesh.star0;
	/*
	const int dimoff[3] = {0,1,1+mesh.n_edges};
	const double bval[3] = {eps_z, 1, 0};
	const complex_t basecoeff[6] = {
		// to higher, to lower
		complex_t(0, 1), complex_t(0, 0), // E->H, none
		complex_t(0,-1), complex_t(0, 1), // H->E, H->divH
		complex_t(0, 0), complex_t(0,-1)  // none, divH->H
	};
	
	// For each matrix index
	for(int ind = 0; ind < nel; ++ind){
		NEWROW();
		const int eldim = ind2el[2*ind+0]; // element dimension
		const int elind = ind2el[2*ind+1]; // element index
		// diagonal
		SETCOL(i, j, ind); SETVAL(-assembly_data.shift*bval[eldim]);
		// Exterior derivative to higher dimension
		if(eldim < 2){
			const int ld = mesh.ld[eldim];
			const unsigned char *d = mesh.d[eldim];
			// For each higher dimension simplex, ...
			for(int el = 0; el < mesh.drows[eldim]; ++c){
				// Check to see if it is incident on this one
				for(int c = 0; c < mesh.dcols[eldim]; ++c){
					// Skip if not
					if(elind != (d[el+c*ld] & LibUMesh_d_ind_MASK)){ continue; }
					int mii[2] = {i,j};
					complex_t coeff = basecoeff[2*eldim+0];
					if(d[el+c*ld] & LibUMesh_d_sign_MASK){ coeff = -coeff; }
					// Fixup each dimension
					for(int d = 0; d < 2; ++d){
						if(0 != (d[el+c*ld] & (LibUMesh_d_uoff_MASK << d))){
							if(0 == mii[d]){
								mii[d] = res[d]-1;
								coeff /= assembly_data.Bloch[d];
							}else{
								mii[d]--;
							}
						}
					}
					SETCOL(mii[0], mii[1], dimoff[eldim+1]+el); SETVAL(coeff);
					// cannot break here, since a single simplex can be incident multiple times
				}
			}
		}
		// Adjoint exterior derivative to lower dimension
		if(eldim > 0){
			const int ld = mesh.ld[eldim-1];
			const unsigned char *d = mesh.d[eldim-1];
			for(int c = 0; c < mesh.dcols[eldim-1]; ++c){
				int mii[2] = {i,j};
				complex_t coeff = basecoeff[2*eldim+1];
				if(d[elind+c*ld] & LibUMesh_d_sign_MASK){ coeff = -coeff; }
				// Fixup each dimension
				for(int d = 0; d < 2; ++d){
					if(0 != (d[elind+c*ld] & (LibUMesh_d_uoff_MASK << d))){
						if(mii[d]+1 == res[d]){
							mii[d] = 0;
							coeff *= assembly_data.Bloch[d];
						}else{
							mii[d]++;
						}
					}
				}
				const int el = d[elind+c*ld] & LibUMesh_d_ind_MASK;
				SETCOL(mii[0], mii[1], dimoff[eldim-1]+el); SETVAL(coeff);
			}
		}
	}
	
	int matbits = matind[q];
	while(matbits & 0xF){
		int m = (matbits & 0xF)-1;
		if(m > 0){
			m--;
			const Material& curmat = material[m];
			int col = nel;
			for(size_t pi = 0; pi < curmat.poles.size(); ++pi){
				const LorentzPole &pole = material[m].poles[pi];
				// P = -i w0 eps V
				NEWROW();
				{
					SETCOL(i, j, col); SETVAL(-assembly_data.shift*eps_z);
					SETCOL(i, j, col+1); SETVAL(complex_t(0.,-pole.omega_0) * eps_z);
				}
				// V = i w0 eps P - i wp eps E
				NEWROW();
				{
					SETCOL(i, j, col+1); SETVAL(-assembly_data.shift*eps_z);
					SETCOL(i, j, col); SETVAL(complex_t(0., pole.omega_0) * eps_z);
					SETCOL(i, j, EZ_OFF); SETVAL(complex_t(0.,-pole.omega_p) * eps_z);
				}
				col += 2;
			}
		}
		matbits >>= 4;
	}
	*/
	
	// Hy = idr[0] * (Ez[i+1,j,k] - Ez[i,j,k])
	//    + idr[1] * (divH[i,j,k] - divH[i,j-1,k])
	NEWROW();
	{
		// diagonal
		SETCOL(i, j, HV_OFF); SETVAL(-assembly_data.shift/mesh.star1[1]);
		// Ez coupling
		if(i+1 == res[0]){
			SETCOL(0, j, EZ_OFF); SETVAL(Im1*assembly_data.Bloch[0]);
		}else{
			SETCOL(i+1, j, EZ_OFF); SETVAL(Im1);
		}
		SETCOL(i, j, EZ_OFF); SETVAL(-Im1);
		// divH coupling
		if(0 == j){
			SETCOL(i, res[1]-1, DIVH_OFF); SETVAL(-Im1/assembly_data.Bloch[1]);
		}else{
			SETCOL(i, j-1, DIVH_OFF); SETVAL(-Im1);
		}
		SETCOL(i, j, DIVH_OFF); SETVAL(Im1);
	}
	// Ez = -idr[1] * (Hx[i,j,k] - Hx[i,j-1,k])
	//    +  idr[0] * (Hy[i,j,k] - Hy[i-1,j,k])
	NEWROW();
	{
		// diagonal
		SETCOL(i, j, EZ_OFF); SETVAL(-assembly_data.shift*eps_z);
		// Hy coupling
		if(0 == i){
			SETCOL(res[0]-1, j, HV_OFF); SETVAL(-Im1/assembly_data.Bloch[0]);
		}else{
			SETCOL(i-1, j, HV_OFF); SETVAL(-Im1);
		}
		SETCOL(i, j, HV_OFF); SETVAL(Im1);
		// Hx coupling
		if(0 == j){
			SETCOL(i, res[1]-1, HU_OFF); SETVAL(Im1/assembly_data.Bloch[1]);
		}else{
			SETCOL(i, j-1, HU_OFF); SETVAL(Im1);
		}
		SETCOL(i, j, HU_OFF); SETVAL(-Im1);
		// Ignore V coupling since we only need the part of the equations
		// to the left of the diagonal.
	}
	// Hx += -idr[1] * (Ez[i,j+1,k] - Ez[i,j,k])
	//    +   idr[0] * (divH[i,j,k] - divH[i-1,j,k])
	NEWROW();
	{
		// diagonal
		SETCOL(i, j, HU_OFF); SETVAL(-assembly_data.shift/mesh.star1[0]);
		// Ez coupling
		if(j+1 == res[1]){
			SETCOL(i, 0, EZ_OFF); SETVAL(-Im1*assembly_data.Bloch[1]);
		}else{
			SETCOL(i, j+1, EZ_OFF); SETVAL(-Im1);
		}
		SETCOL(i, j, EZ_OFF); SETVAL(Im1);
		// divH coupling
		if(0 == i){
			SETCOL(res[0]-1, j, DIVH_OFF); SETVAL(-Im1/assembly_data.Bloch[0]);
		}else{
			SETCOL(i-1, j, DIVH_OFF); SETVAL(-Im1);
		}
		SETCOL(i, j, DIVH_OFF); SETVAL(Im1);
	}
	// divH = idr[0] * (Hx[i+1,j,k] - Hx[i,j,k])
	//      + idr[1] * (Hy[i,j+1,k] - Hy[i,j,k])
	NEWROW();
	{
		// diagonal
		SETCOL(i, j, DIVH_OFF); SETVAL(0.);
		// Hy coupling
		if(j+1 == res[0]){
			SETCOL(i, 0, HV_OFF); SETVAL(Im1*assembly_data.Bloch[1]);
		}else{
			SETCOL(i, j+1, HV_OFF); SETVAL(Im1);
		}
		SETCOL(i, j, HV_OFF); SETVAL(-Im1);
		// Hx coupling
		if(i+1 == res[0]){
			SETCOL(0, j, HU_OFF); SETVAL(Im1*assembly_data.Bloch[0]);
		}else{
			SETCOL(i+1, j, HU_OFF); SETVAL(Im1);
		}
		SETCOL(i, j, HU_OFF); SETVAL(-Im1);
	}
	
	int matbits = matind[q];
	while(matbits & 0xF){
		int m = (matbits & 0xF)-1;
		if(m > 0){
			m--;
			const Material& curmat = material[m];
			int col = 1+mesh.n_edges+mesh.n_faces;
			for(size_t pi = 0; pi < curmat.poles.size(); ++pi){
				const LorentzPole &pole = material[m].poles[pi];
				// P = -i w0 eps V
				NEWROW();
				{
					SETCOL(i, j, col); SETVAL(-assembly_data.shift*eps_z);
					SETCOL(i, j, col+1); SETVAL(complex_t(0.,-pole.omega_0) * eps_z);
				}
				// V = i w0 eps P - i wp eps E
				NEWROW();
				{
					SETCOL(i, j, col+1); SETVAL(-assembly_data.shift*eps_z);
					SETCOL(i, j, col); SETVAL(complex_t(0., pole.omega_0) * eps_z);
					SETCOL(i, j, EZ_OFF); SETVAL(complex_t(0.,-pole.omega_p) * eps_z);
				}
				col += 2;
			}
		}
		matbits >>= 4;
	}
	// finalize with NNZ
	*rowptr = next_col;
	
	assembly_data.p++;
	//SPB_VERB(9, "<GetNextBlockNumeric\n");
	return row_count;
}

double cost_cb(int i, int j, void *data){
	double *cost = (double*)data;
	return cost[i+j*6];
}

int SPB::BandSolver_Ez::MakeMesh(){
	SPB_VERB(9, ">MakeMesh\n");
	SPB_VERB(1, "Preparing mesh\n");
	const size_t Ngrid = res[0] * res[1];

	if(NULL != ind2cell){ SPB_Free(ind2cell, "ind2cell at " __FILE__ ":" STR(__LINE__)); }
	if(NULL != cell2ind){ SPB_Free(cell2ind, "cell2ind at " __FILE__ ":" STR(__LINE__)); }
	if(NULL != matind){ SPB_Free(matind, "matind at " __FILE__ ":" STR(__LINE__)); }
	if(NULL != npoles){ SPB_Free(npoles, "npoles at " __FILE__ ":" STR(__LINE__)); }
	ind2cell = SPB_Alloc(int, 2*Ngrid, "ind2cell at " __FILE__ ":" STR(__LINE__));
	cell2ind = SPB_Alloc(int, Ngrid, "cell2ind at " __FILE__ ":" STR(__LINE__));
	matind = SPB_Alloc(int, Ngrid, "matind at " __FILE__ ":" STR(__LINE__));
	npoles = SPB_Alloc(int, Ngrid, "npoles at " __FILE__ ":" STR(__LINE__));
	
	const double Lr[4] = {
		L.Lr[0], L.Lr[1],
		L.Lr[2], L.Lr[3]
	};
	LibUMesh2_Create(&L.Lr[0], &L.Lr[2], &mesh);
	mesh.star0 /= (double)(res[0]*res[1]);
	
	// Prepare the micro indexing
	{
		nel = 1 + mesh.n_edges + mesh.n_faces;
		double cost[36];
		int eldim[6], eloff[6];
		memset(cost, 0, sizeof(double) * 6*6);
		eldim[0] = 0;
		eloff[0] = 0;
		for(int i = 0; i < mesh.n_edges; ++i){
			cost[0+(1+i)*6] = cost[(1+i)+0*6] = -1;
			eldim[1+i] = 1;
			eloff[1+i] = i;
		}
		for(int f = 0; f < mesh.n_faces; ++f){
			eldim[1+mesh.n_edges+f] = 2;
			eloff[1+mesh.n_edges+f] = f;
			for(int e = 0; e < mesh.dcols[1]; ++e){
				if(0 == (mesh.d21[f+e*2] & (LibUMesh_d_uoff_MASK|LibUMesh_d_voff_MASK))){
					cost[(1+f)+e*6] = -1;
					cost[e+(1+f)*6] = -1;
				}
			}
		}
		int matching[6];
		MinimumWeightPerfectMatching(nel, &cost_cb, matching, cost);
		int j = 0;
		for(int i = 0; i < nel; ++i){
			if(matching[i] < 0){ continue; }
			ind2el[4*j+0] = eldim[i];
			ind2el[4*j+1] = eloff[i];
			ind2el[4*j+2] = eldim[matching[i]];
			ind2el[4*j+3] = eloff[matching[i]];
			matching[matching[i]] = -1;
			j++;
		}
		int dimoff[3] = {0, 1, 1+mesh.n_edges};
		int dimcnt[3] = {0,0,0};
		for(int i = 0; i < nel; ++i){
			int dim = ind2el[2*i+0];
			el2ind[i] = dimoff[dim]+dimcnt[dim];
			dimcnt[dim]++;
		}
	}

	// Prepare the macro indexing
	epsval.clear();
	epsind.resize(Ngrid);
	int max_poles = 0;
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
					num_poles += curmat.poles.size();
					eps_z = curmat.eps_inf.value[8].real();
				}
				npoles[q] = num_poles;
				if(num_poles > max_poles){ max_poles = num_poles; }
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
		for(size_t p = 0; p < Ngrid; ++p){
			int q = ind2cell[2*p+1];
			int i,j; UNIDX(q,i,j);
			ind2cell[2*p+0] = next_index;
			cell2ind[q] = next_index;
			
			// We have Ez, Hx, Hy, divH, plus V,P pairs for each pole
			next_index += 1+mesh.n_edges+mesh.n_faces + 2*npoles[q];
		}
		N = next_index;
	}
	
	assembly_data.max_nnz_per_row = 5+max_poles;
	assembly_data.max_block_size = 1+mesh.n_edges+mesh.n_faces + 2*max_poles;

	ldl.Analyze(*this);
	
	SPB_VERB(9, "<MakeMesh\n");
	return 0;
}

void SPB::BandSolver_Ez::ShiftInv(const complex_t *x, complex_t *y) const{
	RNP::TBLAS::Copy(N, x,1, y,1);
	ldl.Solve(y);
}
void SPB::BandSolver_Ez::Bop(const std::complex<double> *x, std::complex<double> *y) const{
	const int VP_off = 1+mesh.n_edges+mesh.n_faces;
	for(int i = 0; i < res[0]; ++i){
		for(int j = 0; j < res[1]; ++j){
			const int q = IDX(i,j);
			const int row0 = cell2ind[q];
			double eps_z = epsval[epsind[q]]*mesh.star0;
			y[row0+HV_OFF] = x[row0+HV_OFF]/mesh.star1[1];
			y[row0+EZ_OFF] = eps_z*x[row0+EZ_OFF];
			y[row0+HU_OFF] = x[row0+HU_OFF]/mesh.star1[0];
			y[row0+DIVH_OFF] = 0;
			for(int k = 0; k < 2*npoles[q]; ++k){
				y[row0+VP_off+k] = eps_z*x[row0+VP_off+k];
			}
		}
	}
}

void SPB::BandSolver_Ez::PrepareOperator(){
	SPB_Free(cell2ind, "cell2ind at " __FILE__ ":" STR(__LINE__)); cell2ind = NULL;
	SPB_Free(ind2cell, "ind2cell at " __FILE__ ":" STR(__LINE__)); ind2cell = NULL;
	SPB_Free(matind, "matind at " __FILE__ ":" STR(__LINE__)); matind = NULL;
	SPB_Free(npoles, "npoles at " __FILE__ ":" STR(__LINE__)); npoles = NULL;
	MakeMesh();
	
	assembly_data.Bloch[0] = complex_t(cos(k[0]*2*M_PI), sin(k[0]*2*M_PI));
	assembly_data.Bloch[1] = complex_t(cos(k[1]*2*M_PI), sin(k[1]*2*M_PI));
	ldl.Analyze(*this);
}
void SPB::BandSolver_Ez::SetShift(double shift){
	assembly_data.shift = shift;
}
void SPB::BandSolver_Ez::Inertia(int *nlower, int *nupper){
	ldl.Factorize(*this);
	// swapped due to inverse
	ldl.Inertia(nupper, nlower);
}

int SPB::BandSolver_Ez::GetN() const{
	return N;
}
