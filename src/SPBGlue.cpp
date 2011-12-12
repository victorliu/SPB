extern "C" {
#include "SPB.h"
}
#include "SPB.hpp"

// Glue between public C interface and private C++ implementation.

void SPB_BandSolver_Destroy(SPB_BandSolver *S_){
	SPB::BandSolver *S = reinterpret_cast<SPB::BandSolver*>(S_);
	delete S;
}

