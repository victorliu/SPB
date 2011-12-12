#include "SPB.hpp"

SPB::EigenSolver::EigenSolver(const EigenOperator *Op):op(Op),
	data(NULL),
	n_wanted(0),
	target(0),
	tol(1e-7),
	verbosity(9),
	want_interval(0)
{
}

SPB::EigenSolver::~EigenSolver(){
}
