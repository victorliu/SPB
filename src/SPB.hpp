#ifndef _SPB_HPP_INCLUDED_
#define _SPB_HPP_INCLUDED_

#include <string>
#include <vector>
#include <list>
#include <map>
#include <complex>
#include <Sparse.h>
#include "SPB.h"
extern "C" {
#include "ShapeSet.h"
#include "libumesh.h"
}
#include "HermitianMatrixProvider.h"
#include "LDL2.h"

#define SPB_VERB(LVL, ...) do{ \
		if(verbosity >= (LVL)){ \
			fprintf(stdout, __VA_ARGS__); fflush(stdout); \
		} \
	}while(0)

namespace SPB{

typedef std::complex<double> complex_t;

struct ApproximateFrequency{
	double lower, upper;
	int n;
};

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
		const ConstitutiveTensor &eps
	);
	Material(const Material &mat);
};

class Shape{
	int dim;
	virtual void *GetBaseShape() const = 0;
	friend class ShapeSet;
public:
	Shape(int dim):dim(dim){}
	virtual ~Shape(){}
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
	int QueryPt(const double *p, int *tag) const{
		if(2 == dim){
			return ShapeSet2_query_pt(shapeset.ss2, p, NULL, tag);
		}else{
			return ShapeSet3_query_pt(shapeset.ss3, p, NULL, tag);
		}
	}
};


class EigenOperator : public HermitianMatrixProvider{
public:
	// Newly defined by EigenOperator
	virtual void SetShift(double shift) = 0;
	virtual void Inertia(int *nlower, int *nupper) = 0;
	
	virtual void Bop(const complex_t *from, complex_t *to) const = 0;
	virtual void ShiftInv(const complex_t *from, complex_t *to) const = 0;
};

class IntervalEigensolver{
public:
	struct IntervalCount{
		double a, b;
		int n;
	};
	typedef std::list<IntervalCount> interval_list_t;
	typedef HermitianMatrixProvider::complex_t complex_t;

	// set the lower and upper bound of eigenvalues that are sought
	// max_values is the maximum number of eigenvalues to find (set to zero for unlimited)
	// exclude_zero is whether to exclude zero eigenvalues if the range straddles the origin.
	IntervalEigensolver();
	IntervalEigensolver(double lower, double upper, int max_values = 0, bool exclude_zero = true);
	~IntervalEigensolver();
	
	class SolutionHandler{
	public:
		virtual int OnFoundSolution(double value, complex_t *vec) = 0;
	};
	int SetSolutionFunction(SolutionHandler *func);
	class ProgressFunction{
	public:
		virtual int OnProgressUpdate(int nvecs, int nint, double *intbegin, int *intcnt) = 0;
	};
	int SetProgressFunction(ProgressFunction *func);
	
	void SetInterval(double lower, double upper);
	void SetTolerance(double tol);
	
	int SolveCold(EigenOperator *A);
	int SolveWarm(EigenOperator *A, double max_change = 0);
	
	const interval_list_t& GetIntervals() const;
private:
	double range[2];
	int maxvals;
	bool exzero;
	std::list<IntervalCount> ivals;
	
	double tol;
	
	// compute this many eigenvectors at a time
	int block_size;
	
	// maximum number of subintervals of the range to use
	int max_intervals;
	
	SolutionHandler *solfunc;
	ProgressFunction *progfunc;
	
	void AddInterval(double lower, double upper, int count);
	void SearchInterval(EigenOperator *A, double a, double b, int na, int nb);
};

class BandSolver : public EigenOperator{
protected:
	int dim;
	int res[3];
	
	int verbosity;
	
	std::vector<Material> material;
	std::map<std::string,size_t> matmap;
	ShapeSet shapeset;
	
	// Solution
	double k[3];
	double last_k[3]; // in Lk basis
	
	struct{
		complex_t *work;
		size_t n_alloc;
		size_t n_arnoldi;
	} IRA_data;
	size_t n_wanted;
	
	IntervalEigensolver interval_solver;
	LDL2 ldl;
	
	virtual void PrepareOperator(){}
public:
	BandSolver(const Lattice &L);
	virtual ~BandSolver();
	
	void SetResolution(size_t *N);
	void SetVerbosity(int verb){ verbosity = verb; }
	
	void AddMaterial(const Material &mat){ matmap[mat.name] = material.size(); material.push_back(mat); }
	int AddMaterialLorentzPole(const char *name, const LorentzPole &pole);
	int AddShape(const Shape &s, const std::string &matname);
	
	virtual int OutputEpsilon(int *res, const char *filename, const char *format) const = 0;
	
	void ClearSolution();
	void SetK(const double *k);
	int GetApproximateFrequencies(
		double lower, double upper,
		double tol,
		std::list<ApproximateFrequency> &freqs
	);
	int GetBandsNear(
		double target, int n_bands,
		double tol,
		std::vector<double> &freqs,
		complex_t **bands
	);
	int GetPerturbedFrequencies(
		int n_bands,
		std::vector<double> &freqs,
		complex_t **bands
	);
};

class BandSolver_Ez : public BandSolver{
	typedef SPB::complex_t complex_t;
	typedef RNP::Sparse::TCRSMatrix<complex_t> sparse_t;
	Lattice2 L;
	
	size_t N;
	
	// cell2ind is length N.
	// Maps a cell (i,j) with index q=res[1]*i+j to the starting matrix index
	// So cell2ind[q] is the starting matrix index for cell (i,j)
	int *cell2ind;
	
	// ind2cell is length 2*N; pairs of numbers.
	// First in pair maps from permuted cell index to starting matrix index.
	//   This list is in sorted ascending order.
	// Second in pair maps from permuted cell index to cell index
	//   This is the inverse mapping of cell2ind; cell2ind[ind2cell[2*q+0]] == q
	// To map from matrix index back to a cell index, perform a binary search
	//   on ind[2*p+0] to obtain p (round down). Then apply q = ind[2*p+1].
	int *ind2cell;
	
	// matind is length N.
	// Each entry is a bitmap which indicates which materials are present in a cell/
	// matind[q] for cell index q is laid out in blocks of 4 bits (for up to 8 materials
	// in a cell, and up to 15 addressable materials total (zero is regarded as no
	// material), in order from lower order bits to higher order bits until a set
	// of 4 zero bits is encountered.
	int *matind;
	
	// ordering of the simplices within a cell (within a matrix block)
	// first int is the dimension of the simplex
	// second int is the index of the corresponding mesh element.
	int ind2el[2*6];
	int el2ind[6];
	int nel; // number of elements (<= 6)
	
	// npoles is length N.
	// npoles[q] is the number of poles at cell index q.
	int *npoles;

	// vector of all unique epsilon values
	std::vector<int> epsind;
	std::vector<double> epsval;
	
	UMesh2 mesh;
	
	mutable struct{
		complex_t Bloch[2];
		int p;
		int max_block_size;
		int max_nnz_per_row;
		double shift;
	} assembly_data;
	
	void PrepareOperator();
	int MakeMesh();
	int UpdateA(const double k[2]);
	
	void PrintField(const std::complex<double> *x, const char *filename = NULL) const;
protected:
	void Bop(const std::complex<double> *from, std::complex<double> *to) const;
	void ShiftInv(const complex_t *from, complex_t *to) const;
	int GetN() const;
public:
	BandSolver_Ez(double Lr[4]);
	~BandSolver_Ez();
	
	void SetResolution(size_t *N){ res[0] = N[0]; res[1] = N[1]; }
	int OutputEpsilon(int *res, const char *filename, const char *format) const;
	
	void SetShift(double shift);
	void Inertia(int *nlower, int *nupper);
	int GetMaxBlockSize() const;
	int GetMaxNNZPerRow() const;
	int BeginBlockSymbolic() const;
	int GetNextBlockSymbolic(int *rowptr, int *colind) const;
	int BeginBlockNumeric() const;
	int GetNextBlockNumeric(int *rowptr, int *colind, complex_t *value) const;
};

}; // namespace SPB

#endif // _SPB_HPP_INCLUDED_
