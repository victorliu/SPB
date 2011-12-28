#ifndef _INTERVAL_EIGENSOLVER_H_INCLUDED_
#define _INTERVAL_EIGENSOLVER_H_INCLUDED_

#include "HermitianMatrixProvider.h"

class IntervalEigensolver{
public:
	typedef HermitianMatrixProvider::complex_t complex_t;

	// set the lower and upper bound of eigenvalues that are sought
	// max_values is the maximum number of eigenvalues to find (set to zero for unlimited)
	// exclude_zero is whether to exclude zero eigenvalues if the range straddles the origin.
	IntervalEigensolver(double lower, double upper, int max_values = 0, bool exclude_zero = true);
	~IntervalEigensolver();
	
	class SolutionHandler{
	public:
		virtual OnFoundSolution(double value, complex_t *vec) = 0;
	};
	int SetSolutionFunction(const SolutionHandler *func);
	class ProgressFunction{
	public:
		virtual OnProgressUpdate(int nvecs, int nint, double *intbegin, int *intcnt) = 0;
	};
	int SetProgressFunction(const SolutionFunction *func);
	
	int SolveCold(const HermitianMatrixProvider& A);
	int SolveWarm(const HermitianMatrixProvider& A);
private:
	// compute this many eigenvectors at a time
	int block_size = 16;
	
	// maximum number of subintervals of the range to use
	int max_intervals = 8;
	
	SolutionHandler *solfunc;
	ProgressFunction *progfunc;
};

#endif // _INTERVAL_EIGENSOLVER_H_INCLUDED_
