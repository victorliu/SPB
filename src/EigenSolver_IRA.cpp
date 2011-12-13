#include "SPB.hpp"
#include <cstdlib>
#include <IRA.h>
#include <iostream>

typedef SPB::complex_t complex_t;

static void op_(size_t n, const complex_t &shift, const complex_t *x, complex_t *y, void *data){
	const SPB::EigenOperator* op = reinterpret_cast<const SPB::EigenOperator*>(data);
	op->ShiftInv(shift, x, y);
}
static void bv_(size_t n, const complex_t *x, complex_t *y, void *data){
	const SPB::EigenOperator* op = reinterpret_cast<const SPB::EigenOperator*>(data);
	op->Bop(x, y);
}

class SPB::EigenSolver_IRA::Impl{
	size_t n_arnoldi;
	size_t n_alloc;
	friend class SPB::EigenSolver_IRA;
};

SPB::EigenSolver_IRA::EigenSolver_IRA(const EigenOperator *Op):EigenSolver(Op),impl(new Impl()){
	impl->n_arnoldi = 0;
	impl->n_alloc = 0;
}
SPB::EigenSolver_IRA::~EigenSolver_IRA(){
	if(NULL != data){ free(data); }
	delete impl;
}


size_t SPB::EigenSolver_IRA::GetSolutionSize() const{
	return n_wanted;
}
complex_t *SPB::EigenSolver_IRA::GetEigenvalues() const{
	return data;
}
complex_t *SPB::EigenSolver_IRA::GetEigenvectors() const{
	return data + n_wanted;
}
	
void SPB::EigenSolver_IRA::SetNumWanted(size_t k){
	n_wanted = k;
	size_t n = op->GetSize();
	impl->n_arnoldi = 2*k+1;
	size_t n_alloc_new = k + n*impl->n_arnoldi;
	if(impl->n_alloc < n_alloc_new){
		impl->n_alloc = n_alloc_new;
		data = (complex_t*)realloc(data, sizeof(complex_t) * impl->n_alloc);
	}
}

int SPB::EigenSolver_IRA::Solve(){
	size_t n = op->GetSize();
	size_t k = n_wanted;
	SetNumWanted(k);
	SPB::complex_t *w = data;
	SPB::complex_t *v = w+k;
	int nconv = RNP::IRA::ShiftInvert(
		n, target, &op_, &bv_,
		k, impl->n_arnoldi, &RNP::IRA::LargestMagnitude,
		w, v, n,
		NULL,
		NULL,
		(void*)op,
		(void*)op);
	return nconv;
}
