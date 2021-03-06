#include "SPB.hpp"
#include <cstdlib>
#include <IRA.h>

typedef SPB::complex_t complex_t;

SPB::BandSolver::BandSolver(const Lattice &L):dim(L.dim),shapeset(L.dim, L.Lr){
	// set defaults
	res[0] = 16;
	res[1] = 16;
	res[2] = ((2 == dim) ? 1 : 16);
	
	IRA_data.work = NULL;
	IRA_data.n_arnoldi = 0;
	IRA_data.n_alloc = 0;
	
	material.reserve(16);
	SetVerbosity(9);
}
SPB::BandSolver::~BandSolver(){
	ClearSolution();
}

void SPB::BandSolver::ClearSolution(){
}

int SPB::BandSolver::AddShape(const Shape &s, const std::string &matname){
	std::map<std::string,size_t>::const_iterator i = matmap.find(matname);
	if(matmap.end() != i){
		int tag = (int)i->second;
		shapeset.Add(s, tag);
	}else{
		return -2;
	}
	return 0;
}

void SPB::BandSolver::SetResolution(size_t *N){
	res[0] = N[0];
	res[1] = N[1];
	res[2] = ((2 == dim) ? 1 : N[2]);
}

int SPB::BandSolver::AddMaterialLorentzPole(const char *name, const LorentzPole &pole){
	std::map<std::string,size_t>::iterator i = matmap.find(name);
	if(matmap.end() != i){
		material[i->second].poles.push_back(pole);
	}else{
		return -1;
	}
	return 0;
}
/*
void SPB::BandSolver::SetNumBands(size_t k){
	n_wanted = k;
	size_t n = GetProblemSize();
	IRA_data.n_arnoldi = 2*k+1;
	size_t n_alloc_new = k + n*IRA_data.n_arnoldi;
	if(IRA_data.n_alloc < n_alloc_new){
		IRA_data.n_alloc = n_alloc_new;
		IRA_data.work = (complex_t*)realloc(IRA_data.work, sizeof(complex_t) * IRA_data.n_alloc);
	}
}

void SPB::BandSolver::SetTargetFrequencyRange(double lower, double upper){
	target[0] = lower;
	target[1] = upper;
	SetInterval(lower, upper);
}
*/
static void op_(size_t n, const complex_t &shift, const complex_t *x, complex_t *y, void *data){
	const SPB::EigenOperator* op = reinterpret_cast<const SPB::EigenOperator*>(data);
	op->ShiftInv(x, y);
}
static void bv_(size_t n, const complex_t *x, complex_t *y, void *data){
	const SPB::EigenOperator* op = reinterpret_cast<const SPB::EigenOperator*>(data);
	op->Bop(x, y);
}
/*
int SPB::BandSolver::IRASolve(size_t n){
	size_t k = n_wanted;
	SetNumBands(k);
	SPB::complex_t *w = IRA_data.work;
	SPB::complex_t *v = w+k;
	int nconv = RNP::IRA::ShiftInvert(
		n, target[0], &op_, &bv_,
		k, IRA_data.n_arnoldi, &RNP::IRA::LargestMagnitude,
		w, v, n,
		NULL,
		NULL,
		(void*)this,
		(void*)this);
	return nconv;
}

complex_t* SPB::BandSolver::GetFrequencies() const{
	return IRA_data.work;
}
size_t SPB::BandSolver::GetNumSolutions() const{
	return n_wanted;
}*/

void SPB::BandSolver::SetK(const double *newk){
	last_k[0] = k[0];
	k[0] = newk[0];
	last_k[1] = k[1];
	k[1] = newk[1];
	if(dim > 2){
		last_k[2] = k[2];
		k[2] = newk[2];
	}
}
int SPB::BandSolver::GetApproximateFrequencies(
	double lower, double upper,
	double tol,
	std::list<ApproximateFrequency> &freqs
){
	SPB_VERB(1, "Getting approximate frequencies in [%.14g, %.14g]\n", lower, upper);
	ClearSolution();
	PrepareOperator();
	interval_solver.SetInterval(lower, upper);
	interval_solver.SetTolerance(tol);
	interval_solver.SolveCold(this);
	
	freqs.clear();
	SPB::IntervalEigensolver::interval_list_t f = interval_solver.GetIntervals();
	for(SPB::IntervalEigensolver::interval_list_t::const_iterator i = f.begin(); i != f.end(); ++i){
		ApproximateFrequency t;
		t.lower = i->a;
		t.upper = i->b;
		t.n = i->n;
		freqs.push_back(t);
	}
}
