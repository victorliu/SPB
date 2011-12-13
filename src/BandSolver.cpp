#include "SPB.hpp"
#include <cstdlib>

SPB::BandSolver::BandSolver(const Lattice &L):dim(L.dim),shapeset(L.dim, L.Lr){
	// set defaults
	res[0] = 16;
	res[1] = 16;
	res[2] = ((2 == dim) ? 1 : 16);
	
	solver = new SPB::EigenSolver_IRA(this);
	SetTargetFrequency(0);
	SetTolerance(1e-7);
	material.reserve(16);
	SetVerbosity(9);
}
SPB::BandSolver::~BandSolver(){
	ClearSolution();
	delete solver;
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
}

std::vector<SPB::complex_t> SPB::BandSolver::GetFrequencies() const{
	std::vector<SPB::complex_t> ret;
	SPB::complex_t *a = solver->GetEigenvalues();
	for(size_t i = 0; i < solver->GetSolutionSize(); ++i){
		ret.push_back(a[i]);
	}
	return ret;
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
}

