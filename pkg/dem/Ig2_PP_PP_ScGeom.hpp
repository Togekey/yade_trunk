/*CWBoon 2015 */
/* C.W. Boon, G.T. Houlsby, S. Utili (2013).  A new contact detection algorithm for three-dimensional non-spherical particles.  Powder Technology, 248, pp 94-102. */
/* A code for calling MOSEK (ver 6) exists in a previous version of this script. MOSEK can be used to solve the second order cone programme of equations (SOCP), alternativelly to the code used here.*/

#pragma once
#ifdef YADE_POTENTIAL_PARTICLES
// XXX never do #include<Python.h>, see https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/building/include_issues.html
#include <boost/python/detail/wrap_python.hpp>

#include <lib/serialization/Serializable.hpp>
#include <pkg/dem/PotentialParticle.hpp>
#include <pkg/common/Dispatching.hpp>
#include <pkg/common/Sphere.hpp>
#include <Eigen/Core>
#include <stdio.h>

namespace yade { // Cannot have #include directive inside.

class Ig2_PP_PP_ScGeom: public IGeomFunctor {
	public :
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		virtual bool goReverse(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		Real evaluatePP(const shared_ptr<Shape>& cm1, const State& state1, const Vector3r& shift2, const Vector3r newTrial);
		void getPtOnParticle2(const shared_ptr<Shape>& cm1, const State& state1, const Vector3r& shift2, Vector3r previousPt, Vector3r searchDir, Vector3r& newlocalPoint);
		bool customSolve(const shared_ptr<Shape>& cm1, const State& state1, const shared_ptr<Shape>& cm2, const State& state2, const Vector3r& shift2, Vector3r &contactPt, bool warmstart);
		Vector3r getNormal(const shared_ptr<Shape>& cm1, const State& state1, const Vector3r& shift2, const Vector3r newTrial);
		void BrentZeroSurf(const shared_ptr<Shape>& cm1, const State& state1, const Vector3r& shift2, const Vector3r bracketA, const Vector3r bracketB, Vector3r& zero);
		Real getAreaPolygon2(const shared_ptr<Shape>& cm1, const State& state1, const shared_ptr<Shape>& cm2, const State& state2, const Vector3r& shift2, const Vector3r contactPt, const Vector3r contactNormal, int& smaller, Vector3r shearDir, Real& jointLength, const bool twoDimension, Real unitWidth2D, int areaStep);

	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(Ig2_PP_PP_ScGeom,IGeomFunctor,"EXPERIMENTAL. IGeom functor for PotentialParticle - PotentialParticle pair",
			((Real, accuracyTol, pow(10,-7),, "accuracy desired, tolerance criteria for SOCP"))
//			((Real, interactionDetectionFactor, 1.0,,"")) //FIXME: Not used but may be useful in the future for distant interactions
			((Vector3r, twoDdir, Vector3r(0,1,0),, "Direction of 2D"))
			((bool, twoDimension, false,,"Whether the contact is 2-D"))
			((Real, unitWidth2D, 1.0, ,"Unit width in 2D"))
			((bool, calContactArea, true,,"Whether to calculate jointLength for 2-D contacts and contactArea for 2-D and 3-D contacts"))
			((int, areaStep, 5, ,"Angular step (degrees) to calculate :yref:`KnKsPhys.contactArea`. Must be a divisor of 360, e.g. 1,2,3,4,5,6,8,9 and so on, to form a closed loop. Must be smaller than 90 degrees. Smaller angles lead to more accurate calculations but are more expensive"))
			, /* ctor */
		);
	// clang-format on

		FUNCTOR2D(PotentialParticle,PotentialParticle);
		// needed for the dispatcher, even if it is symmetric
		DEFINE_FUNCTOR_ORDER_2D(PotentialParticle,PotentialParticle);
		DECLARE_LOGGER;
};

REGISTER_SERIALIZABLE(Ig2_PP_PP_ScGeom);

} // namespace yade

#ifdef __cplusplus
// note: extern "C" strips out any namespaces.
extern "C" {
#endif

	/* LAPACK LU */
	//int dgesv(int varNo, int varNo2, ::yade::Real *H, int varNo3, int *pivot, ::yade::Real* g, int varNo4, int info){
	extern void dgesv_(const int *N, const int *nrhs, ::yade::Real *Hessian, const int *lda, int *ipiv, ::yade::Real *gradient, const int *ldb, int *info);
	// int ans;
	// dgesv_(&varNo, &varNo2, H, &varNo3, pivot,g, &varNo4, &ans);
	// return ans;
	//}

	/* LAPACK Cholesky */
	extern void dpbsv_(const char *uplo, const int *n, const int *kd, const int *nrhs, ::yade::Real *AB, const int *ldab, ::yade::Real *B, const int *ldb, int *info);

	/* LAPACK QR */
	extern void dgels_(const char *Trans, const int *m, const int *n, const int *nrhs, ::yade::Real *A, const int *lda, ::yade::Real *B, const int *ldb, const ::yade::Real *work, const int *lwork, int *info);


	/*BLAS */
	extern void dgemm_(const char *transA, const char *transB, const int *m, const int *n, const int *k, const ::yade::Real *alpha, ::yade::Real *A, const int *lda, ::yade::Real *B, const int *ldb, const ::yade::Real *beta, ::yade::Real *C, const int *ldc);

	extern void dgemv_(const char *trans, const int *m, const int *n, const ::yade::Real *alpha, ::yade::Real *A, const int *lda, ::yade::Real *x, const int *incx, const ::yade::Real *beta, ::yade::Real *y, const int *incy);

	extern void dcopy_(const int *N, ::yade::Real *x, const int *incx, ::yade::Real *y, const int *incy);

	extern ::yade::Real ddot_(const int *N, ::yade::Real *x, const int *incx, ::yade::Real *y, const int *incy);

	extern void daxpy_(const int *N, const ::yade::Real *da, ::yade::Real *dx, const int *incx, ::yade::Real *dy, const int *incy);

	extern void dscal_(const int *N, const ::yade::Real *alpha, ::yade::Real *x, const int *incx);

	void dsyev_(const char *jobz, const char *uplo, const int *N, ::yade::Real *A, const int *lda, ::yade::Real *W, ::yade::Real *work, int *lwork, int *info);

#ifdef __cplusplus
};
#endif

#endif // YADE_POTENTIAL_PARTICLES

