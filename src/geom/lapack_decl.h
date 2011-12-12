#ifndef _LAPACK_DECL_H_
#define _LAPACK_DECL_H_

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

// AMD Core Math Library requires a significant amount of accomodation.
// None of the LAPACK routines require workspace parameters, and several
// BLAS routines do not enforce const-correctness.

# define FCALL(name) name ## _

double FCALL(ddot)(const int* N, 
                   const double* X, const int* incX, 
                   const double* Y, const int* incY);

double FCALL(dnrm2)(const int* N, 
                    const double* X, const int* incX);

void FCALL(dcopy)(const int* N,
                  const double* X, const int* incX,
                  double* Y, const int* incY);

void FCALL(daxpy)(const int* N,
                  const double* alpha,
                  const double* X, const int* incX,
                  double* Y, const int* incY);

void FCALL(dscal)(const int* N,
                  const double* alpha,
                  double* X, const int* incX);

void FCALL(dgemv)(const char* trans, const int* M, const int* N,
                  const double* alpha,
                  const double* A, const int* lda,
                  const double* X, const int* incX,
                  const double* beta,
                  double* Y, const int* incY);

void FCALL(dgemm)(const char* transA, const char* transB,
                  const int* M, const int* N, const int* K,
                  const double* alpha,
                  const double* A, const int* lda,
                  const double* B, const int* ldb,
                  const double* beta,
                  double* C, const int* ldc);

void FCALL(dgels)(const char *trans, const int *m, const int *n,
                  const int *nrhs, double *a, const int *lda,
                  double *b, const int *ldb,
                  double *work, const int *lwork, int *info);

void FCALL(dgesvd)(const char *jobu, const char *jobvt,
                   const int *m, const int *n,
                   double *a, const int *lda, double *s,
                   double *u, const int *ldu, double *vt, const int *ldvt,
                   double *work, const int *lwork, int *info);

void FCALL(dsysv)(const char *uplo, const int *n, const int *nrhs,
                  double const *a, const int *lda, int *ipiv,
                  double *b, const int *ldb,
                  double *work, const int *lwork, int *info);

void FCALL(dgesv)(const int *n, const int *nrhs,
                  double *a, const int *lda, int *ipiv,
                  double *b, const int *ldb, int *info);

void FCALL(dsyr2k)(const char* uplo, const char* trans,
                   const int* N, const int* K,
                   const double* alpha,
                   const double* A, const int* lda,
                   const double* B, const int* ldb,
                   const double* beta,
                   double* C, const int* ldc);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // _LAPACK_DECL_H_
