#ifndef _LDL2_H_INCLUDED_
#define _LDL2_H_INCLUDED_

#include "HermitianMatrixProvider.h"

// internally all quantities deal with 2x2 complex blocks
// thus all index lists refer to indices of the blocks (length A.GetN()/2)
class LDL2{
public:
	typedef HermitianMatrixProvider::complex_t complex_t;
	
	LDL2();
	~LDL2();
	
	int GetN() const;
	
	int Analyze(const HermitianMatrixProvider& A);
	int Factorize(const HermitianMatrixProvider& A);
	
	int Solve(complex_t *x) const;
	int SolveL(complex_t *x) const;
	int SolveD(complex_t *x) const;
	int SolveLt(complex_t *x) const;
	
	int Inertia(int *pos, int *neg) const;
	
private:
	// n is the size of the matrix (number of block-2x2 rows/cols)
	int n;
	
	// Lp points to beginnings of columns of L (length n+1) in Li and Lx
	// (extra element at end of Lp is number of nonzeros in L)
	int *Lp;
	
	// Parent stores the etree
	// Parent[i] gives the parent node of index i in the tree
	int *Parent;
	
	// Lnz stores the number of nonzeros in each column of L (length n)
	int *Lnz;
	
	// Li stores the row indices of the nonzeros of each column (length nnz)
	// Lx stores the values      of the nonzeros of each column (length nnz)
	//   nnz = Lp[n]
	int *Li;
	complex_t *Lx;
	
	// D stores the inverse of the diagonal factor (length 4*n)
	// D is stored in blocks of 4 corresponding to column-major 2x2 blocks
	complex_t *D;
	
	// If state = 0, everything is NULL
	// If state = 1, we have a symbolic factorization of A
	//               n is set correctly
	//               Lp, Parent, Lnz are allocated
	// If state = 2, we have a numeric factorization of A
	//               Everything should be set.
	// Analyze gets to state 1
	// Factorize gets to state 2
	// Refactorize requires state 1, and reaches state 2
	// Solve requires state 2
	int state;
	
	// Pointer to the previous matrix that was factored or analyzed
	const HermitianMatrixProvider *lastA;
	
	// Internal helper function to deallocate arrays when progressing
	// to a lower state.
	void ReduceState(int state);
};

#endif // _LDL2_H_INCLUDED_
