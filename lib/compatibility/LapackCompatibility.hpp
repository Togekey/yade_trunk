// This file is used during transition from double to Real type.
// TODO: all these functions should be replaced to use Eigen MatrixXr and VectorXr. Then this file should be removed.
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
// not used.
//	extern void dgels_(const char *Trans, const int *m, const int *n, const int *nrhs, double *A, const int *lda, double *B, const int *ldb, const double *work, const int *lwork, int *info);


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

#if (YADE_REAL_BIT != 64)
// when compiling with non-double Real we need converting wrapper functions.
#include <lib/high-precision/Real.hpp>
#include <vector>

/* This is a very simple code that works with Eigen and solves a system of linear equations.
   https://eigen.tuxfamily.org/dox/group__MatrixfreeSolverExample.html
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
int main()
{
	int n = 10;
	Eigen::SparseMatrix<double> A = MatrixXr::Random(n,n).sparseView(0.5,1);
	A = A.transpose()*A;
	VectorXr b(n), x;
	b.setRandom();
	// Solve Ax = b using various iterative solver with matrix-free version:
	{
		Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg;
		cg.compute(A);
		x = cg.solve(b);
		std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
	}
	{
		Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> bicg;
		bicg.compute(A);
		x = bicg.solve(b);
		std::cout << "BiCGSTAB: #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << std::endl;
	}
	{
		Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> gmres;
		gmres.compute(A);
		x = gmres.solve(b);
		std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
	}
	{
		Eigen::DGMRES<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> gmres;
		gmres.compute(A);
		x = gmres.solve(b);
		std::cout << "DGMRES:   #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
	}
	{
		Eigen::MINRES<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> minres;
		minres.compute(A);
		x = minres.solve(b);
		std::cout << "MINRES:   #iterations: " << minres.iterations() << ", estimated error: " << minres.error() << std::endl;
	}
}
*/

inline std::vector<double> toDoubleVec(const ::yade::Real* from, int size) {
	std::vector<double> ret(size,0);
	for(int i = 0 ; i < size ; i++ ) {
		ret[i] = static_cast<double>(from[i]);
	}
	return ret;
}

inline void toRealArrPtr(const std::vector<double>& from, ::yade::Real* to, int size) {
	for(int i = 0 ; i < size ; i++ ) {
		to[i] = static_cast<::yade::Real>(from[i]);
	}
}

// TODO: all these functions should be replaced to use Eigen, MatrixXr and VectorXr.
//       Here are only some workaround/wrappers to use (double precision only) lapack-fortran.
inline void dgesv_(const int *N, const int *nrhs, ::yade::Real *Hessian_Real, const int *lda, int *ipiv, ::yade::Real *gradient_Real, const int *ldb, int *info) {
// DGESV computes the solution to system of linear equations A * X = B for GE matrices
	const int leading_dimension_of_the_array_A = *lda;
	const int leading_dimension_of_the_array_B = *ldb;
	int total_sq_size = leading_dimension_of_the_array_A*leading_dimension_of_the_array_A;

	std::vector<double> Hessian_double  = toDoubleVec(Hessian_Real,total_sq_size);
	std::vector<double> gradient_double = toDoubleVec(gradient_Real,leading_dimension_of_the_array_B);

	dgesv_(N, nrhs, Hessian_double.data() , lda, ipiv, gradient_double.data(), ldb, info);

	toRealArrPtr(Hessian_double,Hessian_Real, total_sq_size);
	toRealArrPtr(gradient_double,gradient_Real, leading_dimension_of_the_array_B);
}

inline void dpbsv_(const char *uplo, const int *n, const int *kd, const int *nrhs, ::yade::Real *AB_Real, const int *ldab, ::yade::Real *B_Real, const int *ldb, int *info) {
// DPBSV computes the solution to system of linear equations A * X = B for OTHER matrices
	const int side = *ldab;
	int total_size = side*side;

	std::vector<double> AB_double = toDoubleVec(AB_Real,total_size);
	std::vector<double> B_double  = toDoubleVec(B_Real,side);
	dpbsv_(uplo, n, kd, nrhs, AB_double.data() , ldab, B_double.data(), ldb, info);
	toRealArrPtr(AB_double,AB_Real,total_size);
	toRealArrPtr( B_double, B_Real,side);
}

inline void dgemm_(const char *transA, const char *transB, const int *m, const int *n, const int *k, const ::yade::Real *alpha_Real, ::yade::Real *A_Real, const int *lda, ::yade::Real *B_Real, const int *ldb, const ::yade::Real *beta_Real, ::yade::Real *C_Real, const int *ldc) {
// DGEMM  performs one of the matrix-matrix operations - FIXME use Eigen instead.
//extern void dgemm_(const char *transA, const char *transB, const int *m, const int *n, const int *k, const double *alpha, double *A, const int *lda, double *B, const int *ldb, const double *beta, double *C, const int *ldc);
	const double alpha_double = static_cast<double>(*alpha_Real);
	const double beta_double = static_cast<double>(*beta_Real);
	std::vector<double> A_double = toDoubleVec(A_Real, (*m) * (*k));
	std::vector<double> B_double = toDoubleVec(B_Real, (*k) * (*n));
	std::vector<double> C_double = toDoubleVec(C_Real, (*n) * (*m));
	dgemm_(transA, transB, m, n, k, &alpha_double, A_double.data(), lda, B_double.data(), ldb, &beta_double, C_double.data(), ldc);
	toRealArrPtr( A_double , A_Real, A_double.size() );
	toRealArrPtr( B_double , B_Real, B_double.size() );
	toRealArrPtr( C_double , C_Real, C_double.size() );
}

inline void dgemv_(const char *trans, const int *m, const int *n, const ::yade::Real *alpha_Real, ::yade::Real *A_Real, const int *lda, ::yade::Real *x_Real, const int *incx, const ::yade::Real *beta_Real, ::yade::Real *y_Real, const int *incy) {
// DGEMV  performs one of the matrix-vector operations - FIXME use Eigen instead.
	if(*incx != 1) throw std::runtime_error("dgemv_ wrapper: incx should be 1");
	if(*incy != 1) throw std::runtime_error("dgemv_ wrapper: incy should be 1");
	const double alpha_double = static_cast<double>(*alpha_Real);
	const double beta_double = static_cast<double>(*beta_Real);
	bool blasN = false;
	if( ((*trans) == 'n') or ((*trans) == 'N')) {
		blasN = true;
	}
	std::vector<double> A_double = toDoubleVec(A_Real, (*m) * (*n) );
	std::vector<double> x_double = toDoubleVec(x_Real, blasN ? (*n) : (*m) );
	std::vector<double> y_double = toDoubleVec(y_Real, blasN ? (*m) : (*n) );
	dgemv_(trans, m, n, &alpha_double, A_double.data(), lda, x_double.data(), incx, &beta_double, y_double.data(), incy);
	toRealArrPtr( A_double , A_Real, A_double.size() );
	toRealArrPtr( x_double , x_Real, x_double.size() );
	toRealArrPtr( y_double , y_Real, y_double.size() );
}

inline void dcopy_(const int *N, ::yade::Real *from_x, const int *incx, ::yade::Real *to_y, const int *incy) {
// DCOPY copies a vector, x, to a vector, y.
	if(*incx != 1) throw std::runtime_error("dcopy_ wrapper: incx should be 1");
	if(*incy != 1) throw std::runtime_error("dcopy_ wrapper: incy should be 1");
	// FIXME: that should be just a simple assignment 
	//   Y=X;
	// with normal VectorXr
	for(int i = 0 ; i < *N ; ++i ) {
		to_y[i] = from_x[i];
	}
}

inline ::yade::VectorXr toVectorXr(const ::yade::Real* from, int size) {
	::yade::VectorXr ret(size);
	for(int i = 0 ; i < size ; i++ ) {
		ret[i] = from[i];
	}
	return ret;
}

inline ::yade::Real ddot_(const int *N, const ::yade::Real *x, const int *incx, const ::yade::Real *y, const int *incy) {
// DDOT forms the dot product of two vectors.
	if(*incx != 1) throw std::runtime_error("ddot_ wrapper: incx should be 1");
	if(*incy != 1) throw std::runtime_error("ddot_ wrapper: incy should be 1");
	::yade::VectorXr xx = toVectorXr(x,*N);
	::yade::VectorXr yy = toVectorXr(y,*N);
	return xx.dot(yy);
}

inline void daxpy_(const int *N, const ::yade::Real *da_Real, ::yade::Real *dx_Real, const int *incx, ::yade::Real *dy_Real, const int *incy) {
	if(*incx != 1) throw std::runtime_error("daxpy_ wrapper: incx should be 1");
	if(*incy != 1) throw std::runtime_error("daxpy_ wrapper: incy should be 1");
	const double da_double = static_cast<double>(*da_Real);
	std::vector<double> dx_double = toDoubleVec(dx_Real, *N );
	std::vector<double> dy_double = toDoubleVec(dy_Real, *N );
	daxpy_(N, &da_double, dx_double.data() , incx, dy_double.data() , incy);
	toRealArrPtr( dx_double , dx_Real, dx_double.size() );
	toRealArrPtr( dy_double , dy_Real, dy_double.size() );
}

inline void dscal_(const int *N, const ::yade::Real *alpha_Real, ::yade::Real *x_Real, const int *incx) {
	if(*incx != 1) throw std::runtime_error("dscal_ wrapper: incx should be 1");
	// FIXME - this should be just an Eigen multiplication by scalar, like:  blasQAs *= scaleFactor; // it unrolls loops, and does everything you want.
	for(int i = 0 ; i < *N ; i++ ) {
		x_Real[i] *= *alpha_Real;
	}
}

inline void dsyev_(const char *jobz, const char *uplo, const int *N, ::yade::Real *A_Real, const int *lda, ::yade::Real *W_Real, ::yade::Real *work_Real, int *lwork, int *info) {
// DSYEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices
	std::vector<double> A_double = toDoubleVec(A_Real, (*lda) * (*N) );
	std::vector<double> W_double = toDoubleVec(W_Real, *N );
	std::vector<double> work_double = toDoubleVec(work_Real, *lwork );
	dsyev_(jobz, uplo, N, A_double.data(), lda, W_double.data(), work_double.data(), lwork, info);
	toRealArrPtr( A_double , A_Real, A_double.size() );
	toRealArrPtr( W_double , W_Real, W_double.size() );
	toRealArrPtr( work_double , work_Real, work_double.size() );

}

#endif

