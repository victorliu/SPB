#ifndef _SPB_HPP_INCLUDED_
#define _SPB_HPP_INCLUDED_

#include <string>
#include <vector>
#include <map>
#include <complex>
#include <Sparse.h>
#include "SPB.h"
extern "C" {
#include "ShapeSet.h"
}
#include "HermitianMatrixProvider.h"
#include "IntervalEigensolver.h"
#include "LDL2.h"

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


class EigenOperator{
public:
	virtual size_t GetSize() const = 0;
	virtual void SetShift(double shift) = 0;
	
	virtual void Bop(const complex_t *from, complex_t *to) const = 0;
	virtual void ShiftInv(const complex_t *from, complex_t *to) const = 0;
};

class BandSolver : public EigenOperator, public IntervalEigensolver{
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
	
	struct{
		complex_t *work;
		size_t n_alloc;
		size_t n_arnoldi;
	} IRA_data;
	size_t n_wanted;
	double approxtol, tol;
	double target[2];
	
	virtual size_t GetProblemSize() const = 0;
	virtual void StructureChanged(){};
public:
	BandSolver(const Lattice &L);
	virtual ~BandSolver();
	
	void SetResolution(size_t *N);
	void SetNumBands(size_t n);
	void SetTargetFrequencyRange(double lower, double upper);
	void SetApproximationTolerance(double tolerance){ approxtol = tolerance; }
	void SetTolerance(double tolerance){ tol = tolerance; }
	void SetVerbosity(int verb){ verbosity = verb; }
	
	void AddMaterial(const Material &mat){ matmap[mat.name] = material.size(); material.push_back(mat); }
	int AddMaterialLorentzPole(const char *name, const LorentzPole &pole);
	int AddShape(const Shape &s, const std::string &matname);
	
	virtual int OutputEpsilon(int *res, const char *filename, const char *format) const = 0;
	
	void ClearSolution();
	virtual int SolveK(const double *k) = 0;
	complex_t* GetFrequencies() const;
	size_t GetNumSolutions() const;
	
	int IRASolve(size_t n);
};

class BandSolver_Ez : public BandSolver, public HermitianMatrixProvider{
	typedef SPB::complex_t complex_t;
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
	
	// npoles is length N.
	// npoles[q] is the number of poles at cell index q.
	int *npoles;

	// vector of all unique epsilon values
	std::vector<int> epsind;
	std::vector<double> epsval;

	LDL2 ldl;
	
	struct{
		int n_edges; // 2 or 3
		int edge[3*2*2];
		// edge:
		//        from to
		//         u v u v
		// u vec [ 0 0 1 0 ]
		// v vec [ 0 0 0 1 ]
		// w vec [ 0 0 1 1 ] w = u+v in this case
		// w could also be [ 1 0 0 1 ], w = from u to v = v-u
		
		double star_mu[3];
		// star_mu is hodge star for each edge
		double star_eps; // same as primal cell area
		int face[5*2];
		// face:
		// for a square lattice:
		//        sign
		//   u0 [  1 ]
		//   v0 [ -1 ]
		//   w0 [  0 ]
		//   u1 [ -1 ]
		//   v1 [  1 ]
		//   second set of 5 are all zero
		// for a u.v < 0 lattice:
		//        sign
		//   u0 [  1 ]
		//   v0 [  0 ]
		//   w0 [ -1 ]
		//   u1 [  0 ]
		//   v1 [  1 ]
		//
		//   u0 [  0 ]
		//   v0 [ -1 ]
		//   w0 [  1 ]
		//   u1 [ -1 ]
		//   v1 [  0 ]
		
	} mesh;
	
	mutable struct{
		complex_t Bloch[2];
		int p;
		int max_block_size;
		int max_nnz_per_row;
		double shift;
	} assembly_data;
	
	// When structure changes, invalidate A, B, ind
	// When K changes, symbolic A still good, need to update it
	// When shift changes, symbolic A still good, numeric needs regenerating
	// If ind is not NULL, we assume we have an A that is symbolically factored, and B is valid
	bool valid_A_numeric;
	
	void StructureChanged();
	int MakeMesh();
	int UpdateA(const double k[2]);
	
	void PrintField(const std::complex<double> *x, const char *filename = NULL) const;
protected:
	size_t GetSize() const;
	void Bop(const std::complex<double> *from, std::complex<double> *to) const;
	void ShiftInv(const complex_t *from, complex_t *to) const;
	int GetN() const;
	size_t GetProblemSize() const;
public:
	BandSolver_Ez(double Lr[4]);
	~BandSolver_Ez();
	
	void SetResolution(size_t *N){ res[0] = N[0]; res[1] = N[1]; }
	int OutputEpsilon(int *res, const char *filename, const char *format) const;
	int SolveK(const double *k);
		
	void SetShift(double shift);
	int GetMaxBlockSize() const;
	int GetMaxNNZPerRow() const;
	int BeginBlockSymbolic() const;
	int GetNextBlockSymbolic(int *rowptr, int *colind) const;
	int BeginBlockNumeric() const;
	int GetNextBlockNumeric(int *rowptr, int *colind, complex_t *value) const;
};

}; // namespace SPB

#endif // _SPB_HPP_INCLUDED_
