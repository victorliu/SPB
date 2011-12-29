#include "LDL2.h"
#include <Sparse.h>
#include <cstdio>

typedef HermitianMatrixProvider::complex_t complex_t;
typedef RNP::Sparse::TCRSMatrix<complex_t> sparse_t;

class MatrixLoader : public HermitianMatrixProvider{
	sparse_t *A;
	mutable int iiter;
public:
	MatrixLoader(const char *filename):iiter(0){
		sparse_t::entry_map_t Amap;
		int n, i, j;
		double zr, zi;
		FILE *fp = fopen(filename, "rb");
		fscanf(fp, "%d", &n);
		while(1 == fscanf(fp, "%d", &i)){
			fscanf(fp, "%d", &j);
			fscanf(fp, "%lf", &zr);
			fscanf(fp, "%lf", &zi);
			printf("Adding %d,%d = %g,%g\n", i,j,zr,zi);
			Amap[sparse_t::index_t(i,j)] = complex_t(zr,zi);
		}
		fclose(fp);
		
		A = new sparse_t(n, n, Amap);
	}
	~MatrixLoader(){
		delete A;
	}
	
	int GetN() const{ return A->n; }
	int GetMaxBlockSize() const{ return 2; }
	int GetMaxNNZPerRow() const{ return GetN(); }
	
	int BeginBlockSymbolic() const{ iiter = 0; return 0; }
	int GetNextBlockSymbolic(int *rowptr, int *colind) const{
		if(iiter >= GetN()){ return 0; }
		rowptr[0] = A->rowptr[iiter];
		rowptr[1] = A->rowptr[iiter+1];
		rowptr[2] = A->rowptr[iiter+2];
		for(int r = 0; r < 2; ++r){
			for(int p = rowptr[r]; p < rowptr[r+1]; ++p){
				*colind = A->colind[p];
				colind++;
			}
		}
		rowptr[1] -= rowptr[0];
		rowptr[2] -= rowptr[0];
		rowptr[0] = 0;
		iiter += 2;
		return 2;
	}
	int BeginBlockNumeric() const{ iiter = 0; return 0; }
	int GetNextBlockNumeric(int *rowptr, int *colind, complex_t *value) const{
		if(iiter >= GetN()){ return 0; }
		rowptr[0] = A->rowptr[iiter];
		rowptr[1] = A->rowptr[iiter+1];
		rowptr[2] = A->rowptr[iiter+2];
		for(int r = 0; r < 2; ++r){
			for(int p = rowptr[r]; p < rowptr[r+1]; ++p){
				*colind = A->colind[p];
				colind++;
				*value = A->values[p];
				value++;
			}
		}
		rowptr[1] -= rowptr[0];
		rowptr[2] -= rowptr[0];
		rowptr[0] = 0;
		iiter += 2;
		return 2;
	}
	
	void Mult(const complex_t *x, complex_t *y) const{
		RNP::Sparse::MultMV<'N'>(*A, x, y);
	}
};

int main(int argc, char *argv[]){
	LDL2 ldl;
	char *filename = "A.txt";
	if(argc > 1){
		filename = argv[1];
	}
	MatrixLoader A(filename);
	const int N = A.GetN();
	complex_t *x, *y, *z;
	
	x = new complex_t[N];
	y = new complex_t[N];
	z = new complex_t[N];
	
	for(int i = 0; i < N; ++i){
		double f = ((double)i+0.5) / (double)N;
		y[i] = x[i] = complex_t(f, 1.-f);
	}
	
	ldl.Factorize(A);
	ldl.Solve(y);
	
	A.Mult(y, z);
	
	// y = inv(A)*x
	// z = A*inv(A)*x = x
	
	double norm = 0;
	for(int i = 0; i < N; ++i){
		z[i] -= x[i];
		double diff = std::abs(z[i]);
		if(diff > norm){ norm = diff; }
	}
	
	printf("rnorm = %.14g\n", norm);
	
	delete [] z;
	delete [] y;
	delete [] x;
	return 0;
}
