#include <cstdlib>
extern "C" {
#include "SPB.h"
}
#include "SPB.hpp"

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
		SPB_BandSolver* ret = (SPB_BandSolver*)malloc(sizeof(SPB_BandSolver));
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
	free(S);
}


int SPB_BandSolver_GetDimension(SPB_BandSolver *S){
	if(NULL == S){ return 0; }
	return S->dim;
}
char SPB_BandSolver_GetPolarization(SPB_BandSolver *S){
	if(NULL == S){ return 0; }
	return S->pol;
}
