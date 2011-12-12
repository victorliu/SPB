#ifndef _SPB_HPP_INCLUDED_
#define _SPB_HPP_INCLUDED_

#include <string>
#include <vector>
#include <map>
#include <complex>
#include "SPB.h"
extern "C" {
#include "ShapeSet.h"
}

#define SPB_VERB(LVL, STR, ...) do{ \
		if(verbosity >= LVL){ \
			fprintf(stdout, STR, __VA_ARGS__); fflush(stdout); \
		} \
	}while(0)

namespace SPB{

typedef std::complex<double> complex_t;


class Lattice{
protected:
	Lattice(int dim):dim(dim){}
public:
	int dim;
	double Lr[9];
	double Lk[9];
	virtual ~Lattice(){}
	virtual double CharacteristicKLength() const = 0;
};

class Lattice2 : public Lattice{
public:
	Lattice2(double L[4]);
	double CharacteristicKLength() const;
};
class Lattice3 : public Lattice{
public:
	Lattice3(double L[9]);
	double CharacteristicKLength() const;
};


struct LorentzPole{
	double omega_0, Gamma, omega_p;
	LorentzPole();
	LorentzPole(const LorentzPole &lp);
};

struct ConstitutiveTensor{
	enum{
		SCALAR,
		DIAGONAL,
		TENSOR
	} type;
	complex_t value[9];
	ConstitutiveTensor();
	ConstitutiveTensor(const ConstitutiveTensor &et);
};


struct Material{
	std::string name;
	ConstitutiveTensor eps_inf;
	std::vector<LorentzPole> poles;
	Material(
		const std::string &name,
		const ConstitutiveTensor &eps,
		const std::vector<LorentzPole> &poles
	);
	Material(const Material &mat);
};

class Shape{
	int dim;
	virtual void *GetBaseShape() const = 0;
	friend class ShapeSet;
public:
	Shape(int dim):dim(dim){}
	virtual ~Shape();
	virtual bool Contains(const double *p) = 0;
};

class Shape2 : public Shape{
	shape2 s;
	void *GetBaseShape() const{ return (void*)&s; }
public:
	Shape2(const shape2 &s2):Shape(2),s(s2){}
	~Shape2(){}
	bool Contains(const double *p){
		return (bool)shape2_contains(&s, p);
	}
};

class Shape3 : public Shape{
	shape3 s;
	void *GetBaseShape() const{ return (void*)&s; }
public:
	Shape3(const shape3 &s3):Shape(3),s(s3){}
	~Shape3(){}
	bool Contains(const double *p){
		return (bool)shape3_contains(&s, p);
	}
};

class ShapeSet{
	int dim;
	union{
		ShapeSet2 ss2;
		ShapeSet3 ss3;
	} shapeset;
public:
	ShapeSet(int dim, const double *Lr):dim(dim){
		if(2 == dim){
			shapeset.ss2 = ShapeSet2_new(Lr);
		}else{
			shapeset.ss3 = ShapeSet3_new(Lr);
		}
	}
	~ShapeSet(){
		if(2 == dim){
			ShapeSet2_destroy(shapeset.ss2);
		}else{
			ShapeSet3_destroy(shapeset.ss3);
		}
	}
	int Add(const Shape &s, int tag){
		if(2 == s.dim){
			return ShapeSet2_add(shapeset.ss2, (shape2*)s.GetBaseShape(), tag);
		}else{
			return ShapeSet3_add(shapeset.ss3, (shape3*)s.GetBaseShape(), tag);
		}
	}
	int QueryPt(const double *p, int *tag){
		if(2 == dim){
			return ShapeSet2_query_pt(shapeset.ss2, p, NULL, tag);
		}else{
			return ShapeSet3_query_pt(shapeset.ss3, p, NULL, tag);
		}
	}
};


class EigenOperator{
public:
	virtual size_t GetSize() const = 0;
	
	// Applies the matrices A and B
	virtual void Aop(const complex_t *from, complex_t *to) const = 0;
	virtual void Bop(const complex_t *from, complex_t *to) const = 0;
	
	////// JDQZ specific
	// Preconditioner for A-target*B (not currently used)
	virtual void Precond(const std::complex<double> &alpha, const std::complex<double> &beta, const std::complex<double> *from, std::complex<double> *to) const = 0;
	// Generates a random vector in the range of inv(A)*B (not currently used)
	virtual void Randvec(std::complex<double> *x) const = 0;
	// Performs orthogonalization against nullspace of B. (not currently used)
	virtual void Orth(complex_t *x) const = 0;
	
	////// IRA specific
	virtual void ShiftInv(const complex_t &shift, const complex_t *from, complex_t *to) const = 0;
	// ShiftOp is the inverse of ShiftInv (so, ShiftOp is the forward shifted operator A-shift*B
	// We include this so that implementers are (hopefully) forced to implement the forward FFT-based
	// operation to compare against manually applying Aop and shift*Bop.
//	virtual void ShiftOp(const complex_t &shift, const complex_t *from, complex_t *to) const = 0;
};

class EigenSolver{
protected:
	complex_t *data;
	const EigenOperator *op;
	
	size_t n_wanted;
	complex_t target;
	int want_interval;
	complex_t target2; // upper end of interval if want_interval == 1
	double tol;
	int verbosity;
public:
	EigenSolver(const EigenOperator *Op);
	virtual ~EigenSolver();
	virtual int Solve() = 0;
	virtual size_t GetSolutionSize() const = 0;
	virtual complex_t *GetEigenvalues() const = 0;
	virtual complex_t *GetEigenvectors() const = 0;
	
	virtual void SetNumWanted(size_t nwanted){ n_wanted = nwanted; }
	void SetTarget(const complex_t &targ){ target = targ; }
	void SetTolerance(double tolerance){ tol = tolerance; }
	void SetVerbosity(int verb){ verbosity = verb; }
};

class EigenSolver_JDQZ : public EigenSolver{
	class Impl;
	Impl *impl;
public:
	EigenSolver_JDQZ(const EigenOperator *Op);
	~EigenSolver_JDQZ();
	int Solve();
	size_t GetSolutionSize() const;
	complex_t *GetEigenvalues() const;
	complex_t *GetEigenvectors() const;
	
	void SetNumWanted(size_t nwanted);
};

class EigenSolver_IRA : public EigenSolver{
	class Impl;
	Impl *impl;
public:
	EigenSolver_IRA(const EigenOperator *Op);
	~EigenSolver_IRA();
	int Solve();
	size_t GetSolutionSize() const;
	complex_t *GetEigenvalues() const;
	complex_t *GetEigenvectors() const;
	
	void SetNumWanted(size_t nwanted);
};

class BandSolver : public EigenOperator{
protected:
	int dim;
	int res[3];
	
	int verbosity;
	// method options
	int use_direct_solver; // set to 1 if the shift-invert operator should be directly factored
	int force_hermitian; // set to 1 if we should assume the problem is hermitian
	int force_invsym; // set to 1 if we shoud assume the structure is inversion symmetric
	
	std::vector<Material> material;
	std::map<std::string,size_t> matmap;
	ShapeSet shapeset;
	
	// Solution
	double last_k[3]; // in Lk basis
	EigenSolver *solver;
public:
	BandSolver(const Lattice &L);
	virtual ~BandSolver();
	
	void SetResolution(size_t *N);
	void SetNumBands(size_t n){ solver->SetNumWanted(n); }
	void SetTargetFrequency(const complex_t &freq){ solver->SetTarget(freq); }
	void SetTolerance(double tol){ solver->SetTolerance(tol); }
	void SetVerbosity(int verb){ solver->SetVerbosity(verbosity = verb); }
	
	void AddMaterial(const Material &mat){ matmap[mat.name] = material.size(); material.push_back(mat); }
	int AddShape(const Shape &s, const std::string &matname);
	
	void ClearSolution();
	virtual int SolveK(const double *k) = 0;
	std::vector<complex_t> GetFrequencies() const;
};

class BandSolver_Ez : public BandSolver{
	Lattice2 L;
	class Impl;
	Impl *impl;
	
	void PrintField(const std::complex<double> *x, const char *filename = NULL) const;
protected:
	size_t GetSize() const;
	void Aop(const std::complex<double> *from, std::complex<double> *to) const;
	void Bop(const std::complex<double> *from, std::complex<double> *to) const;
	void Precond(const std::complex<double> &alpha, const std::complex<double> &beta, const std::complex<double> *from, std::complex<double> *to) const;
	void Randvec(complex_t *to) const;
	void Orth(complex_t *x) const;
	void ShiftInv(const complex_t &shift, const complex_t *from, complex_t *to) const;
public:
	BandSolver_Ez(double Lr[4]);
	~BandSolver_Ez();
	
	void SetResolution(size_t *N){ res[0] = N[0]; res[1] = N[1]; }
	
	int SolveK(const double *k);
	
	void Op(size_t n, const complex_t &shift, const complex_t *from, complex_t *to) const;
	void OpForw(size_t n, const complex_t &shift, const complex_t *from, complex_t *to) const;
};

class BandSolver_Hz : public BandSolver{
	Lattice2 L;
	class Impl;
	Impl *impl;
	
	void PrintField(const std::complex<double> *x, const char *filename = NULL) const;
protected:
	size_t GetSize() const;
	void Aop(const std::complex<double> *from, std::complex<double> *to) const;
	void Bop(const std::complex<double> *from, std::complex<double> *to) const;
	void Precond(const std::complex<double> &alpha, const std::complex<double> &beta, const std::complex<double> *from, std::complex<double> *to) const;
	void Randvec(std::complex<double> *x) const;
	void Orth(complex_t *x) const;
	void ShiftInv(const complex_t &shift, const complex_t *from, complex_t *to) const;
public:
	BandSolver_Hz(double Lr[4]);
	~BandSolver_Hz();
	
	void SetResolution(size_t *N){ res[0] = N[0]; res[1] = N[1]; }
	
	int SolveK(const double *k);
};

}; // namespace SPB

#endif // _SPB_HPP_INCLUDED_
