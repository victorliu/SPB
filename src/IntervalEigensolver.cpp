#include "IntervalEigensolver.h"

IntervalEigensolver::IntervalEigensolver(){
	tol = 1e-4;
	block_size = 16;
	max_intervals = 8;
}

IntervalEigensolver::IntervalEigensolver(double lower, double upper, int max_values, bool exclude_zero):
maxvals(max_values),
exzero(exclude_zero)
{
	range[0] = lower;
	range[1] = upper;
}

IntervalEigensolver::~IntervalEigensolver(){
}

int IntervalEigensolver::SetSolutionFunction(SolutionHandler *func){
	solfunc = func;
	return 0;
}
int IntervalEigensolver::SetProgressFunction(ProgressFunction *func){
	progfunc = func;
	return 0;
}

void IntervalEigensolver::SetInterval(double lower, double upper){
	range[0] = lower;
	range[1] = upper;
}

int IntervalEigensolver::SolveCold(const HermitianMatrixProvider& A){
	return 0;
}
	