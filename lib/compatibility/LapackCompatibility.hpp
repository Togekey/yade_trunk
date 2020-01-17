#pragma once

#ifdef __cplusplus
// note: extern "C" strips out any namespaces.
extern "C" {
#endif

/* LAPACK LU */
	//int dgesv(int varNo, int varNo2, double *H, int varNo3, int *pivot, double* g, int varNo4, int info){
	extern void dgesv_(const int *N, const int *nrhs, double *Hessian, const int *lda, int *ipiv, double *gradient, const int *ldb, int *info);
	// int ans;
	// dgesv_(&varNo, &varNo2, H, &varNo3, pivot,g, &varNo4, &ans);
	// return ans;
	//}

/* LAPACK Cholesky */
	extern void dpbsv_(const char *uplo, const int *n, const int *kd, const int *nrhs, double *AB, const int *ldab, double *B, const int *ldb, int *info);

/* LAPACK QR */
	extern void dgels_(const char *Trans, const int *m, const int *n, const int *nrhs, double *A, const int *lda, double *B, const int *ldb, const double *work, const int *lwork, int *info);


/*BLAS */
	extern void dgemm_(const char *transA, const char *transB, const int *m, const int *n, const int *k, const double *alpha, double *A, const int *lda, double *B, const int *ldb, const double *beta, double *C, const int *ldc);

	extern void dgemv_(const char *trans, const int *m, const int *n, const double *alpha, double *A, const int *lda, double *x, const int *incx, const double *beta, double *y, const int *incy);

	extern void dcopy_(const int *N, double *x, const int *incx, double *y, const int *incy);

	extern double ddot_(const int *N, double *x, const int *incx, double *y, const int *incy);

	extern void daxpy_(const int *N, const double *da, double *dx, const int *incx, double *dy, const int *incy);

	extern void dscal_(const int *N, const double *alpha, double *x, const int *incx);

	void dsyev_(const char *jobz, const char *uplo, const int *N, double *A, const int *lda, double *W, double *work, int *lwork, int *info);

#ifdef __cplusplus
};
#endif

#ifdef __cplusplus
extern "C" {
#endif
	void dgesv_(const int *N, const int *nrhs, double *Hessian, const int *lda, int *ipiv, double *gradient, const int *ldb, int *info);
	void dsyev_(const char *jobz, const char *uplo, const int *N, double *A, const int *lda, double *W, double *work, int *lwork, int *info);
#ifdef __cplusplus
};
#endif

#if (YADE_REAL_BIT != 64)
// when compiling with non-double Real we need converting wrapper functions.
#include <lib/high-precision/Real.hpp>

#endif

