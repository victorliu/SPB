#include "LDL2.h"
#include <cstdlib>
extern "C"{
#include "mem.h"
}
//#include <iostream>

LDL2::LDL2():
	Lp(NULL),Parent(NULL),Lnz(NULL),
	Li(NULL),Lx(NULL),D(NULL),
	state(0),n(0),
	lastA(NULL)
{
}

LDL2::~LDL2(){
	ReduceState(0);
}
void LDL2::ReduceState(int targstate){
	if(targstate < 2){
		if(NULL != Li){ SPB_Free(Li, "Li at " __FILE__ ":" STR(__LINE__)); Li = NULL; }
		if(NULL != Lx){ SPB_Free(Lx, "Lx at " __FILE__ ":" STR(__LINE__)); Lx = NULL; }
		if(NULL != D){ SPB_Free(D, "D at " __FILE__ ":" STR(__LINE__)); D = NULL; }
		if(targstate < 1){
			if(NULL != Lp){ SPB_Free(Lp, "Lp at " __FILE__ ":" STR(__LINE__)); Lp = NULL; }
			if(NULL != Parent){ SPB_Free(Parent, "Parent at " __FILE__ ":" STR(__LINE__)); Parent = NULL; }
			if(NULL != Lnz){ SPB_Free(Lnz, "Lnz at " __FILE__ ":" STR(__LINE__)); Lnz = NULL; }
			lastA = NULL;
		}
	}
}

int LDL2::GetN() const{ return n; }

/*
a + bi   c + di
[ a -b ] [ c ]
[ b  a ] [ d ]
*/

// Each 2x2 complex block can be expanded as 4 2x2 real blocks:
// [ z1 z3 ] [ a ]    [ z1r -z3i z3r -z1i ] [ ar ]
// [ z2 z4 ] [ b ] -> [ z2i  z4r z4i  z2r ] [ bi ]
//                    [ z2r -z4i z4r -z2i ] [ br ]
//                    [ z1i  z3r z3i  z1r ] [ ai ]

int LDL2::Analyze(const HermitianMatrixProvider& A){
	ReduceState(0);
	
	n = A.GetN()/2;
	const int nbmax2 = A.GetMaxBlockSize();
	const int nbmax = nbmax2/2;
	
	Lp = SPB_Alloc(int, n+1, "Lp at " __FILE__ ":" STR(__LINE__));
	Parent = SPB_Alloc(int, n, "Parent at " __FILE__ ":" STR(__LINE__));
	Lnz = SPB_Alloc(int, n, "Lnz at " __FILE__ ":" STR(__LINE__));
	
	int *Flag = SPB_Alloc(int, n, "Flag at " __FILE__ ":" STR(__LINE__));
	int *rowptr = SPB_Alloc(int, nbmax2+1, "rowptr at " __FILE__ ":" STR(__LINE__));
	int *colind = SPB_Alloc(int, A.GetMaxNNZPerRow()*nbmax2, "colind at " __FILE__ ":" STR(__LINE__));
	
	A.BeginBlockSymbolic();
	
	int k = 0; // column index of original matrix, also column index of block 2x2 real
	int nr;
	while(nr = A.GetNextBlockSymbolic(rowptr, colind)){
		if(nr <= 0){
			// error
			break;
		}else if(1 == nr % 2){
			// bad; need to have even number of block rows
			break;
		}else{
			nr /= 2;
		}
		for(int k1 = 0; k1 < nr; ++k1, ++k){
			Parent[k] = -1;
			Flag[k] = k;
			Lnz[k] = 0;
			for(int b = 0; b < 2; ++b){
				for(int p = rowptr[2*k1+b]; p < rowptr[2*k1+b+1]; ++p){
					//int i = A.PermInv(colind[p])/2;
					int i = colind[p]/2;
					if(i < k){
						// follow path from i to root of etree, stop at flagged node
						for(; Flag[i] != k; i = Parent[i]){
							// find parent of i if not yet determined
							if(Parent[i] < 0){ Parent[i] = k; }
							Lnz[i]++;				// L (k,i) is nonzero
							Flag[i] = k;			// mark i as visited
						}
					}
				}
			}
		}
	}
	
	/* construct Lp index array from Lnz column counts */
	Lp[0] = 0;
	for(k = 0; k < n; k++){
		Lp[k+1] = Lp[k] + Lnz[k];
	}
	
	SPB_Free(colind, "colind at " __FILE__ ":" STR(__LINE__));
	SPB_Free(rowptr, "rowptr at " __FILE__ ":" STR(__LINE__));
	SPB_Free(Flag, "Flag at " __FILE__ ":" STR(__LINE__));
	
	state = 1;
	lastA = &A;
	return 0;
}

#define C2ZERO(A) do{ \
		(A)[0] = 0.; (A)[1] = 0.; \
		(A)[2] = 0.; (A)[3] = 0.; \
	}while(0)
#define C2SET(A,B) do{ \
		(A)[0] = (B)[0]; (A)[1] = (B)[1]; \
		(A)[2] = (B)[2]; (A)[3] = (B)[3]; \
	}while(0)
#define C2SETH(A,B) do{ \
		(A)[0] = std::conj((B)[0]); (A)[1] = std::conj((B)[2]); \
		(A)[2] = std::conj((B)[1]); (A)[3] = std::conj((B)[3]); \
	}while(0)

// A -= B*C
#define C2SUBPROD(A,B,C) do{ \
		(A)[0] -= (B)[0] * (C)[0]; \
		(A)[0] -= (B)[2] * (C)[1]; \
		(A)[1] -= (B)[1] * (C)[0]; \
		(A)[1] -= (B)[3] * (C)[1]; \
		(A)[2] -= (B)[0] * (C)[2]; \
		(A)[2] -= (B)[2] * (C)[3]; \
		(A)[3] -= (B)[1] * (C)[2]; \
		(A)[3] -= (B)[3] * (C)[3]; \
	}while(0)
// C = A*B
#define C2MUL(A,B,C) do{ \
		(C)[0] = (A)[0] * (B)[0]; \
		(C)[0]+= (A)[2] * (B)[1]; \
		(C)[1] = (A)[1] * (B)[0]; \
		(C)[1]+= (A)[3] * (B)[1]; \
		(C)[2] = (A)[0] * (B)[2]; \
		(C)[2]+= (A)[2] * (B)[3]; \
		(C)[3] = (A)[1] * (B)[2]; \
		(C)[3]+= (A)[3] * (B)[3]; \
	}while(0)
#define C2INV(A) do{ \
		complex_t idet(1./((A)[0]*(A)[3]-(A)[1]*(A)[2])); \
		complex_t t((A)[0]); \
		(A)[0] = (A)[3] * idet; \
		(A)[3] = t * idet; \
		(A)[1] = -(A)[1] * idet; \
		(A)[2] = -(A)[2] * idet; \
	}while(0)

int LDL2::Factorize(const HermitianMatrixProvider& A){
	ReduceState(1);
	if(state < 1 || &A != lastA){
		ReduceState(0);
		Analyze(A);
	}
	
	const int nbmax2 = A.GetMaxBlockSize();
	const int nbmax = nbmax2/2;

	Li = SPB_Alloc(int, Lp[n], "Li at " __FILE__ ":" STR(__LINE__));
	Lx = SPB_Alloc(complex_t, 4*Lp[n], "Lx at " __FILE__ ":" STR(__LINE__));
	D  = SPB_Alloc(complex_t, 4*n, "D at " __FILE__ ":" STR(__LINE__));
	
	complex_t *Y = SPB_Alloc(complex_t, 4*n, "Y at " __FILE__ ":" STR(__LINE__));
	int *Pattern = SPB_Alloc(int, n, "Pattern at " __FILE__ ":" STR(__LINE__));
	int *Flag = SPB_Alloc(int, n, "Flag at " __FILE__ ":" STR(__LINE__));
	int *rowptr = SPB_Alloc(int, nbmax2+1, "rowptr at " __FILE__ ":" STR(__LINE__));
	int *colind = SPB_Alloc(int, A.GetMaxNNZPerRow()*nbmax2, "colind at " __FILE__ ":" STR(__LINE__));
	complex_t *colval = SPB_Alloc(complex_t, A.GetMaxNNZPerRow()*nbmax2, "colval at " __FILE__ ":" STR(__LINE__));
	
	A.BeginBlockNumeric();
	
	int k = 0; // column index of original matrix, also column index of block 2x2 real
	int nr;
	int len;
	while(nr = A.GetNextBlockNumeric(rowptr, colind, colval)){
		nr /= 2;
		for(int k1 = 0; k1 < nr; ++k1, ++k){
			C2ZERO(&Y[4*k]); // Y(0:k) is now all zero
			int top = n;     // stack for pattern is empty
			Flag[k] = k;     // mark node k as visited
			Lnz[k] = 0;      // count of nonzeros in column k of L
			for(int b = 0; b < 2; ++b){
				for(int p = rowptr[2*k1+b]; p < rowptr[2*k1+b+1]; ++p){
					//int i_ = A.PermInv(colind[p]);
					int i_ = colind[p];
					int i = i_/2;
					int ir= i_%2;
					if(i <= k){
						// scatter A(i,k) into Y (sum duplicates)
						Y[4*i+2*b+ir] += std::conj(colval[p]);
						for(len = 0; Flag[i] != k; i = Parent[i]){
							Pattern[len++] = i;	// L(k,i) is nonzero
							Flag[i] = k;        // mark i as visited
						}
						while(len > 0) Pattern[--top] = Pattern[--len];
					}
				}
			}
			// compute numerical values kth row of L (a sparse triangular solve)
			C2SET(&D[4*k], &Y[4*k]); // get D(k,k) and clear Y(k)
			C2ZERO(&Y[4*k]);
			for(; top < n; top++){
				int i = Pattern[top]; // Pattern[top:n-1] is pattern of L(:,k)
				complex_t yi[4];
				complex_t Di[4];
				complex_t l_ki[4];
				int p2 = Lp[i] + Lnz[i];
				C2SET(yi, &Y[4*i]); // get and clear Y(i)
				C2ZERO(&Y[4*i]);
				int p;
				for(p = Lp[i]; p < p2; p++){
					C2SUBPROD(&Y[4*Li[p]], &Lx[4*p], yi);
				}
				C2SET(Di, &D[4*i]);
				C2INV(Di);
				C2MUL(Di, yi, l_ki); // the nonzero entry L(k,i)

				C2SETH(&Lx[4*p], l_ki);
				C2SUBPROD(&D[4*k], &Lx[4*p], yi);
				Li[p] = k; // store L(k,i) in column form of L
				Lnz[i]++;  //increment count of nonzeros in col i
			}
		}
	}
	/*
	for(k = 0; k < n; k++){
		//std::cout << "[ " << D[4*k+0] << ", " << D[4*k+2] << std::endl;
		//std::cout << "  " << D[4*k+1] << ", " << D[4*k+3] << " ]" << std::endl;
		C2INV(&D[4*k]);
	}*/
	
	SPB_Free(colval, "colval at " __FILE__ ":" STR(__LINE__));
	SPB_Free(colind, "colind at " __FILE__ ":" STR(__LINE__));
	SPB_Free(rowptr, "rowptr at " __FILE__ ":" STR(__LINE__));
	SPB_Free(Flag, "Flag at " __FILE__ ":" STR(__LINE__));
	SPB_Free(Pattern, "Pattern at " __FILE__ ":" STR(__LINE__));
	SPB_Free(Y, "Y at " __FILE__ ":" STR(__LINE__));
	
	state = 2;
	lastA = &A;
	return 0;
}
int LDL2::Solve(complex_t *x) const {
	int ret;
	ret = SolveL(x);
	if(0 != ret){ return ret; }
	ret = SolveD(x);
	if(0 != ret){ return ret; }
	ret = SolveLt(x);
	if(0 != ret){ return ret; }
	return 0;
}

int LDL2::SolveL(complex_t *x) const{
	if(state < 2){ return -1; }

	for(int j = 0; j < n; j++){
		int p2 = Lp[j+1];
		for(int p = Lp[j]; p < p2; p++){
			//X[Li[p]] -= Lx[p] * X[j];
			x[2*Li[p]+0] -= Lx[4*p+0] * x[2*j+0];
			x[2*Li[p]+0] -= Lx[4*p+2] * x[2*j+1];
			x[2*Li[p]+1] -= Lx[4*p+1] * x[2*j+0];
			x[2*Li[p]+1] -= Lx[4*p+3] * x[2*j+1];
		}
	}
	return 0;
}
int LDL2::SolveD(complex_t *x) const{
	if(state < 2){ return -1; }

	/*
	for(int j = 0; j < n; j++){
		//X[j] *= D[j]; // D is actually inv(D)
		complex_t s, t;
		s = D[4*j+0] * x[2*j+0];
		s+= D[4*j+2] * x[2*j+1];
		t = D[4*j+1] * x[2*j+0];
		t+= D[4*j+3] * x[2*j+1];
		x[2*j+0] = s;
		x[2*j+1] = t;
	}*/
	for(int j = 0; j < n; j++){
		//X[j] *= D[j]; // D is actually inv(D)
		complex_t det(D[4*j+0]*D[4*j+3]-D[4*j+1]*D[4*j+2]);
		complex_t s((D[4*j+3] * x[2*j+0] - D[4*j+2] * x[2*j+1]) / det);
		complex_t t((D[4*j+0] * x[2*j+1] - D[4*j+1] * x[2*j+0]) / det);
		x[2*j+0] = s;
		x[2*j+1] = t;
	}
	return 0;
}
int LDL2::SolveLt(complex_t *x) const{
	if(state < 2){ return -1; }

	for(int j = n-1; j >= 0; j--){
		int p2 = Lp[j+1];
		for(int p = Lp[j]; p < p2; p++){
			//X[j] -= Lx[p] * X[Li[p]];
			x[2*j+0] -= std::conj(Lx[4*p+0]) * x[2*Li[p]+0];
			x[2*j+0] -= std::conj(Lx[4*p+1]) * x[2*Li[p]+1];
			x[2*j+1] -= std::conj(Lx[4*p+2]) * x[2*Li[p]+0];
			x[2*j+1] -= std::conj(Lx[4*p+3]) * x[2*Li[p]+1];
		}
	}
	return 0;
}
int LDL2::Inertia(int *pos, int *neg) const{
	if(NULL == pos){ return -1; }
	if(NULL == neg){ return -2; }
	if(state < 2){ return -3; }
	
	*pos = 0;
	*neg = 0;
	for(int j = 0; j < n; j++){
		double ah = 0.5 * D[4*j+0].real();
		double ch = 0.5 * D[4*j+3].real();
		double ac = ah+ch;
		double bmag = std::abs(D[4*j+1]);
		double root = hypot(ah-ch, bmag);
		if(ac + root > 0){
			(*pos)++;
		}else if(ac + root < 0){
			(*neg)++;
		}
		if(ac - root > 0){
			(*pos)++;
		}else if(ac - root < 0){
			(*neg)++;
		}
	}
	return 0;
}
