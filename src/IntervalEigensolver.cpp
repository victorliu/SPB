#include "SPB.hpp"
#include <stack>
#include <list>

#include <iostream>
#include <cstdlib>

SPB::IntervalEigensolver::IntervalEigensolver(){
	tol = 1e-2;
	block_size = 16;
	max_intervals = 8;
}

SPB::IntervalEigensolver::IntervalEigensolver(double lower, double upper, int max_values, bool exclude_zero):
maxvals(max_values),
exzero(exclude_zero)
{
	range[0] = lower;
	range[1] = upper;
}

SPB::IntervalEigensolver::~IntervalEigensolver(){
}

int SPB::IntervalEigensolver::SetSolutionFunction(SolutionHandler *func){
	solfunc = func;
	return 0;
}
int SPB::IntervalEigensolver::SetProgressFunction(ProgressFunction *func){
	progfunc = func;
	return 0;
}

void SPB::IntervalEigensolver::SetInterval(double lower, double upper){
	range[0] = lower;
	range[1] = upper;
}
void SPB::IntervalEigensolver::SetTolerance(double tolerance){
	tol = tolerance;
}


struct ival{
	double a, b;
	int na, nb;
};
void SPB::IntervalEigensolver::SearchInterval(SPB::EigenOperator *A, double a, double b, int na, int nb){
	int dummy;
	ival newival;
	newival.a = a;
	newival.b = b;
	newival.na = na;
	newival.nb = nb;
	
	std::stack<ival> istack;
	istack.push(newival);
	while(!istack.empty()){
		ival top = istack.top();
		istack.pop();
//std::cout << "a,b=" << top.a << "," << top.b << ", na,nb=" << top.na << "," << top.nb << std::endl;
		if(top.nb <= top.na){ continue; }
		if(top.b - top.a < tol*(range[1]-range[0])){
			IntervalCount c;
			c.a = top.a;
			c.b = top.b;
			c.n = top.nb-top.na;
//std::cout << "  [" << c.a << "," << c.b << "]:" << c.n << std::endl;
			ivals.push_back(c);
			continue;
		}
		double m = 0.5*(top.a + top.b);
		int nm;
		A->SetShift(m);
		A->Inertia(&nm, &dummy);
		newival.a = m; newival.b = top.b;
		newival.na = nm; newival.nb = top.nb;
		istack.push(newival);
		newival.a = top.a; newival.b = m;
		newival.na = top.na; newival.nb = nm;
		istack.push(newival);
	}
}

/*
struct IntervalCountSorter{
	bool operator()(
		const IntervalEigensolver::IntervalCount &a,
		const IntervalEigensolver::IntervalCount &b
	) const{
		return a.a < b.a;
	}
};*/
int SPB::IntervalEigensolver::SolveCold(SPB::EigenOperator *A){
	int na, nb, dummy;
	/*
	for(double s = 0.01; s <= 10.0; s += 0.1){
		A->SetShift(s);
		A->Inertia(&na, &dummy);
		std::cerr << "s = " << s/(2*M_PI) << ", lower = " << na << ", upper = " << dummy << std::endl;
	}
	exit(0);
*/
	ivals.clear();
	A->SetShift(range[0]);
	A->Inertia(&na, &dummy);
std::cerr << "lower = " << na << ", upper = " << dummy << std::endl;
	A->SetShift(range[1]);
	A->Inertia(&nb, &dummy);
	
	SearchInterval(A, range[0], range[1], na, nb);
	
	//std::sort(ivals.begin(), ivals.end(), IntervalCountSorter());
	return 0;
}

int SPB::IntervalEigensolver::SolveWarm(SPB::EigenOperator *A, double max_change){
	int na, nb, dummy;
	double a, b;
	// For each interval, expand by maxchange and re-analyze
	// Start by making sure that the regions outside are still as they should be
	int nlower;
	int alower = range[0];
	A->SetShift(alower);
	A->Inertia(&nlower, &dummy);
	const std::list<IntervalCount> icopy(ivals);
	ivals.clear();
	for(std::list<IntervalCount>::const_iterator i = icopy.begin(); i != icopy.end(); ++i){
		a = i->a;
		b = i->b;
		const int n = i->n;
		// expand current interval by max_change
		a -= max_change;
		b += max_change;
		if(a < alower){
			a = alower;
			na = nlower;
			A->SetShift(b);
			A->Inertia(&nb, &dummy);
		}else{
			A->SetShift(a);
			A->Inertia(&na, &dummy);
			A->SetShift(b);
			A->Inertia(&nb, &dummy);
			if(na > nlower){
				SearchInterval(A, alower, a, nlower, na);
			}
		}
		SearchInterval(A, a, b, na, nb);
		nlower = nb;
		alower = b;
	}
	if(b < range[1]){
		int nupper;
		A->SetShift(range[1]);
		A->Inertia(&nupper, &dummy);
		if(nb > nlower){
			SearchInterval(A, b, range[1], nb, nupper);
		}
	}
	return 0;
}

const SPB::IntervalEigensolver::interval_list_t& SPB::IntervalEigensolver::GetIntervals() const{
	return ivals;
}

