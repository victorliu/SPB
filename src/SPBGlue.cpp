#include <cstdlib>
#include <cstring>
extern "C" {
#include "SPB.h"
}
#include "SPB.hpp"
#include "util.h"
extern "C"{
#include "mem.h"
}

// Glue between public C interface and private C++ implementation.

struct tag_SPB_BandSolver{
	int dim;
	char pol;
	SPB::BandSolver *S;
};

SPB_BandSolver* SPB_BandSolver_New(int dim, char pol, double *Lr){
	SPB::BandSolver *S = NULL;
	if(2 == dim){
		if('E' == pol || 'e' == pol){
			S = new SPB::BandSolver_Ez(Lr);
		}else if('H' == pol || 'h' == pol){
		}
	}else if(3 == dim){
	}
	if(NULL != S){
		SPB_BandSolver* ret = SPB_Alloc(SPB_BandSolver, 1, "ret at " __FILE__ ":" STR(__LINE__));
		if(NULL != ret){
			ret->dim = dim;
			ret->pol = pol;
			ret->S = S;
			return ret;
		}
	}
	return NULL;
}

void SPB_BandSolver_Destroy(SPB_BandSolver *S){
	if(NULL == S){ return; }
	delete S->S;
	SPB_Free(S, "S at " __FILE__ ":" STR(__LINE__));
}


int SPB_BandSolver_GetDimension(const SPB_BandSolver *S){
	if(NULL == S){ return 0; }
	return S->dim;
}
char SPB_BandSolver_GetPolarization(const SPB_BandSolver *S){
	if(NULL == S){ return 0; }
	return S->pol;
}





int SPB_BandSolver_AddMaterial(SPB_BandSolver *S, const char *name, const SPB_ConstitutiveTensor *eps){
	if(NULL == S){ return -1; }
	if(NULL == name){ return -2; }
	if(NULL == eps){ return -3; }
	SPB::ConstitutiveTensor eps_inf;
	switch(eps->type){
	case SPB_ConstitutiveTensor::SPB_ConstitutiveTensor_SCALAR:
		eps_inf.type = SPB::ConstitutiveTensor::SCALAR;
		break;
	case SPB_ConstitutiveTensor::SPB_ConstitutiveTensor_DIAGONAL:
		eps_inf.type = SPB::ConstitutiveTensor::DIAGONAL;
		break;
	case SPB_ConstitutiveTensor::SPB_ConstitutiveTensor_TENSOR:
		eps_inf.type = SPB::ConstitutiveTensor::TENSOR;
		break;
	}
	memcpy(eps_inf.value, eps->value, sizeof(SPB::complex_t) * 9);
	S->S->AddMaterial(SPB::Material(name, eps_inf));
	return 0;
}
int SPB_BandSolver_SetMaterial(SPB_BandSolver *S, const char *name, const SPB_ConstitutiveTensor *eps){
	if(NULL == S){ return -1; }
	if(NULL == name){ return -2; }
	if(NULL == eps){ return -3; }
	return -1;
}
int SPB_BandSolver_Material_AddLorentzPole(SPB_BandSolver *S, const char *name, const SPB_LorentzPole *pole){
	if(NULL == S){ return -1; }
	if(NULL == name){ return -2; }
	if(NULL == pole){ return -3; }
	SPB::LorentzPole p;
	p.omega_0 = pole->omega_0;
	p.omega_p = pole->omega_p;
	p.Gamma = pole->Gamma;
	S->S->AddMaterialLorentzPole(name, p);
	
	return 0;
}
int SPB_BandSolver_RemoveMaterial(SPB_BandSolver *S, const char *name){
	if(NULL == S){ return -1; }
	if(NULL == name){ return -2; }
	return -1;
}

int SPB_BandSolver_AddRectangle(SPB_BandSolver *S,
	const char *material,
	double center[2],
	double halfwidth[2],
	double angle
){
	shape2 r;
	double ca, sa, idet;
	if(NULL == S){ return -1; }
	if(2 != S->dim){ return -1; }
	if(NULL == material){ return -2; }
	if(NULL == center){ return -2; }
	if(NULL == halfwidth){ return -2; }
	
	ca = cos(angle);
	sa = sin(angle);
	
	r.type = SHAPE2_QUAD;
	r.o[0] = center[0];
	r.o[1] = center[1];
	r.quad.A[0] = ca*halfwidth[0];
	r.quad.A[1] = sa*halfwidth[0];
	r.quad.A[2] = -sa*halfwidth[1];
	r.quad.A[3] = ca*halfwidth[1];
	idet = 1./(halfwidth[0]*halfwidth[1]);
	r.quad.B[0] = idet * ca*halfwidth[1];
	r.quad.B[1] = idet *-sa*halfwidth[0];
	r.quad.B[2] = idet * sa*halfwidth[1];
	r.quad.B[3] = idet * ca*halfwidth[0];
	S->S->AddShape(SPB::Shape2(r), material);
	
	return -1;
}

int SPB_BandSolver_OutputEpsilon(const SPB_BandSolver *S,
	int *res,
	const char *filename,
	const char *format
){
	if(NULL == S){ return -1; }
	if(NULL == res){ return -2; }
	S->S->OutputEpsilon(res, filename, format);
	return -1;
}

int SPB_BandSolver_SetResolution(SPB_BandSolver *S, int *res){
	int i;
	size_t n[3];
	if(NULL == S){ return -1; }
	if(NULL == res){ return -2; }
	for(i = 0; i < S->dim; ++i){
		n[i] = res[i];
		if(n[i] < 1){ return -2; }
	}
	S->S->SetResolution(n);
	return 0;
}
int SPB_BandSolver_SetVerbosity(SPB_BandSolver *S, int v){
	if(NULL == S){ return -1; }
	S->S->SetVerbosity(v);
	return 0;
}
int SPB_BandSolver_SetK(SPB_BandSolver *S, double *k){
	if(NULL == S){ return -1; }
	if(NULL == k){ return -2; }
	S->S->SetK(k);
	return 0;
}

int SPB_BandSolver_GetApproximateFrequencies(
	SPB_BandSolver *S,
	double lower, double upper,
	double tol,
	int *n, SPB_ApproximateFrequency **lst
){
	if(NULL == S){ return -1; }
	if(lower <= 0){ return -2; }
	if(upper <= lower){ return -3; }
	if(tol <= 0){ return -4; }
	if(NULL == n){ return -5; }
	if(NULL == lst){ return -6; }
	
	std::list<SPB::ApproximateFrequency> freqs;
	S->S->GetApproximateFrequencies(lower, upper, tol, freqs);
	
	*lst = NULL;
	for(std::list<SPB::ApproximateFrequency>::const_iterator i = freqs.begin(); i != freqs.end(); ++i){
		*lst = SPB_Alloc(SPB_ApproximateFrequency, 1, "lst at " __FILE__ ":" STR(__LINE__));
		(*lst)->lower = i->lower;
		(*lst)->upper = i->upper;
		(*lst)->n = i->n;
		(*lst)->next = NULL;
		lst = &((*lst)->next);
	}
	return 0;
}
