/*CWBoon 2016 */
/* Please cite: */
/* CW Boon, GT Houlsby, S Utili (2013).  A new algorithm for contact detection between convex polygonal and polyhedral particles in the discrete element method.  Computers and Geotechnics 44, 73-82. */
/* The numerical library is changed from CPLEX to CLP because subscription to the academic initiative is required to use CPLEX for free */

#ifdef YADE_POTENTIAL_BLOCKS

#pragma once
// XXX never do #include<Python.h>, see https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/building/include_issues.html
#include <boost/python/detail/wrap_python.hpp>

#include <lib/serialization/Serializable.hpp>
#include <pkg/dem/PotentialBlock.hpp>
#include <pkg/common/Dispatching.hpp>
#include <pkg/common/Sphere.hpp>
#include <Eigen/Core>
#include <stdio.h>

#include <ClpSimplex.hpp>
#include <CoinHelperFunctions.hpp>
#include <CoinTime.hpp>
#include <CoinBuild.hpp>
#include <CoinModel.hpp>

#include <iomanip>
#include <cassert>

namespace yade { // Cannot have #include directive inside.

class Ig2_PB_PB_ScGeom: public IGeomFunctor
{
	//protected:

	public :
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		virtual bool goReverse(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		Real getSignedArea(const Vector3r pt1,const Vector3r pt2, const Vector3r pt3);
		Real evaluatePB(const shared_ptr<Shape>& cm1, const State& state1, const Vector3r& shift2, const Vector3r newTrial);
		void getPtOnParticle2(const shared_ptr<Shape>& cm1, const State& state1, const Vector3r& shift2, Vector3r previousPt, Vector3r searchDir, Vector3r& newlocalPoint);
//		void getPtOnParticleArea(const shared_ptr<Shape>& cm1, const State& state1, const Vector3r& shift2, Vector3r previousPt, Vector3r normal, Vector3r& newlocalPoint);
		bool getPtOnParticleAreaNormal(const shared_ptr<Shape>& cm1, const State& state1, const Vector3r& shift2, const Vector3r previousPt, const Vector3r prevDir, const int prevNo, Vector3r& newlocalPoint, Vector3r& newNormal, int& newNo);
		Real getDet(const MatrixXr A);
		bool customSolve(const shared_ptr<Shape>& cm1, const State& state1, const shared_ptr<Shape>& cm2, const State& state2, const Vector3r& shift2, Vector3r &contactPt, bool warmstart);
		Real evaluatePhys(const shared_ptr<Shape>& cm1,  const State& state1, const Vector3r& shift2, const Vector3r newTrial, Real& phi_b, Real& phi_r, /* Real& JRC, Real& JSC, */ Real& cohesion, /* Real& ks, Real& kn, */ Real& tension, /* Real &lambda0, Real &heatCapacity, Real &hwater, */ bool &intactRock, int &activePlanesNo, int &jointType);
		Vector3r getNormal(const shared_ptr<Shape>& cm1, const State& state1, const Vector3r& shift2, const Vector3r newTrial, const bool twoDimension);
		void BrentZeroSurf(const shared_ptr<Shape>& cm1, const State& state1, const Vector3r& shift2, const Vector3r bracketA, const Vector3r bracketB, Vector3r& zero);
		bool startingPointFeasibilityCLP(const shared_ptr<Shape>& cm1, const State& state1, const shared_ptr<Shape>& cm2, const State& state2, const Vector3r& shift2, Vector3r &contactPoint/*, bool &convergeFeasibility*/);
		bool customSolveAnalyticCentre(const shared_ptr<Shape>& cm1, const State& state1, const shared_ptr<Shape>& cm2, const State& state2, const Vector3r& shift2, Vector3r& contactPt);
		Real getAreaPolygon2(const shared_ptr<Shape>& cm1, const State& state1, const shared_ptr<Shape>& cm2, const State& state2, const Vector3r& shift2, const Vector3r contactPt, const Vector3r contactNormal, int& smaller, Vector3r shearDir, Real& jointLength, const bool twoDimension, Real unitWidth2D);


	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(Ig2_PB_PB_ScGeom,IGeomFunctor,"PB",
		((Real, accuracyTol, pow(10,-7),, "accuracy desired, tolerance criteria for SOCP"))
//		((Real, stepAngle, pow(10,-2),, ""))
//		((Real,interactionDetectionFactor,1.0,,"")) //FIXME: Not used but may be useful in the future for distant interactions
		((Vector3r, twoDdir, Vector3r(0,1,0),, "Direction of 2D"))
		((bool, twoDimension, false,,"Whether the contact is 2-D"))
		((Real, unitWidth2D, 1.0,,"Unit width in 2D"))
		((bool, calContactArea, true,,"Whether to calculate jointLength for 2-D contacts and contactArea for 2-D and 3-D contacts"))
		, /* ctor */
	);
	// clang-format on

	FUNCTOR2D(PotentialBlock,PotentialBlock);
	// needed for the dispatcher, even if it is symmetric
	DEFINE_FUNCTOR_ORDER_2D(PotentialBlock,PotentialBlock);
	DECLARE_LOGGER;

};

REGISTER_SERIALIZABLE(Ig2_PB_PB_ScGeom);


#ifdef __cplusplus
extern "C" {
#endif

/* LAPACK LU */
	//int dgesv(int varNo, int varNo2, Real *H, int varNo3, int *pivot, Real* g, int varNo4, int info){
	 extern void dgesv_(const int *N, const int *nrhs, Real *Hessian, const int *lda, int *ipiv, Real *gradient, const int *ldb, int *info);
	// int ans;
	// dgesv_(&varNo, &varNo2, H, &varNo3, pivot,g, &varNo4, &ans);
	// return ans;
	//}

/* LAPACK Cholesky */
	extern void dpbsv_(const char *uplo, const int *n, const int *kd, const int *nrhs, Real *AB, const int *ldab, Real *B, const int *ldb, int *info);

/* LAPACK QR */
	extern void dgels_(const char *Trans, const int *m, const int *n, const int *nrhs, Real *A, const int *lda, Real *B, const int *ldb, const Real *work, const int *lwork, int *info);


/*BLAS */
	extern void dgemm_(const char *transA, const char *transB, const int *m, const int *n, const int *k, const Real *alpha, Real *A, const int *lda, Real *B, const int *ldb, const Real *beta, Real *C, const int *ldc);

	extern void dgemv_(const char *trans, const int *m, const int *n, const Real *alpha, Real *A, const int *lda, Real *x, const int *incx, const Real *beta, Real *y, const int *incy);

	extern void dcopy_(const int *N, Real *x, const int *incx, Real *y, const int *incy);

	extern Real ddot_(const int *N, Real *x, const int *incx, Real *y, const int *incy);

	extern void daxpy_(const int *N, const Real *da, Real *dx, const int *incx, Real *dy, const int *incy);

	extern void dscal_(const int *N, const Real *alpha, Real *x, const int *incx);


	void dsyev_(const char *jobz, const char *uplo, const int *N, Real *A, const int *lda, Real *W, Real *work, int *lwork, int *info);


#ifdef __cplusplus
};
#endif


} // namespace yade

#endif // YADE_POTENTIAL_BLOCKS

