#ifndef _HERMITIAN_MATRIX_PROVIDER_H_INCLUDED
#define _HERMITIAN_MATRIX_PROVIDER_H_INCLUDED

#include <complex>

// Interface of a class that is able to provide blocks of rows of a
// Hermitian matrix. A user of this interface retrieves the blocks
// by first calling BeginBlock___() to allow the implementer to
// initialize the iteration. Then, GetNextBlock___() is called to
// retreive the blocks until there are none left.
class HermitianMatrixProvider{
public:
	typedef std::complex<double> complex_t;
	virtual int GetN() const = 0;
	virtual int GetMaxBlockSize() const = 0;
	virtual int GetMaxNNZPerRow() const = 0;
	
	// By default, there is no permutation
	virtual int Perm(int i) const{ return i; };
	virtual int PermInv(int i) const{ return i; };
	
	// Returns number of rows, assumed to always be even
	// Assumes that diagonal 2x2 blocks are non-singular
	virtual int BeginBlockSymbolic() const = 0;
	virtual int GetNextBlockSymbolic(int *rowptr, int *colind) const = 0;
	virtual int BeginBlockNumeric() const = 0;
	virtual int GetNextBlockNumeric(int *rowptr, int *colind, complex_t *value) const = 0;
};

#endif // _HERMITIAN_MATRIX_PROVIDER_H_INCLUDED
