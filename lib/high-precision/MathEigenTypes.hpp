/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef REAL_EXPOSED_TYPES_HPP
#define REAL_EXPOSED_TYPES_HPP

#include <Eigen/Core>
#include <complex>
#include <unsupported/Eigen/AlignedVector3>

#ifndef YADE_REAL_MATH_NAMESPACE
#error "This file cannot be included alone, include Real.hpp instead"
#endif

namespace yade {

// integral type for indices, to avoid compiler warnings with int
typedef Eigen::Matrix<int, 1, 1>::Index Index;
typedef Eigen::Matrix<int, 2, 1>        Vector2i;
typedef Eigen::Matrix<int, 3, 1>        Vector3i;
typedef Eigen::Matrix<int, 6, 1>        Vector6i;

/*************************************************************************/
/*************************        Real          **************************/
/*************************************************************************/

/* exposed types */
typedef Eigen::Matrix<Real, 2, 1>                           Vector2r;
#ifdef EIGEN_DONT_ALIGN
typedef Eigen::Matrix<Real, 3, 1>                           Vector3r;
#else
typedef Eigen::AlignedVector3<Real>                         Vector3r;
#endif
typedef Eigen::AlignedVector3<Real>                         Vector3ra;
typedef Eigen::Matrix<Real, 4, 1>                           Vector4r;
typedef Eigen::Matrix<Real, 6, 1>                           Vector6r;
typedef Eigen::Matrix<Real, 3, 3>                           Matrix3r;
typedef Eigen::Matrix<Real, 6, 6>                           Matrix6r;
typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> MatrixXr;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 1>              VectorXr;
typedef Eigen::Quaternion<Real>                             Quaternionr;
typedef Eigen::AngleAxis<Real>                              AngleAxisr;
typedef Eigen::AlignedBox<Real, 3>                          AlignedBox3r;
typedef Eigen::AlignedBox<Real, 2>                          AlignedBox2r;

/*************************************************************************/
/*************************       Complex        **************************/
/*************************************************************************/

typedef Eigen::Matrix<Complex, 2, 1>                           Vector2cr;
typedef Eigen::Matrix<Complex, 3, 1>                           Vector3cr;
typedef Eigen::Matrix<Complex, 6, 1>                           Vector6cr;
typedef Eigen::Matrix<Complex, Eigen::Dynamic, 1>              VectorXcr;
typedef Eigen::Matrix<Complex, 3, 3>                           Matrix3cr;
typedef Eigen::Matrix<Complex, 6, 6>                           Matrix6cr;
typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> MatrixXcr;

/*************************************************************************/
/*************************   for external use   **************************/
/*************************************************************************/

// This is for external applications, shouldn't be normally used.
// Use `using namespace ::yade::AllMathTypes;` only inside a .cpp file! Otherwise the types will leak outside which will cause compilation errors due to ambiguity.
namespace MathEigenTypes {
	using ::yade::Real;

	using ::yade::Complex;

	using ::yade::Index;
	using ::yade::Vector2i;
	using ::yade::Vector2r;
	using ::yade::Vector3i;
	using ::yade::Vector3r;
	using ::yade::Vector3ra;
	using ::yade::Vector4r;
	using ::yade::Vector6i;
	using ::yade::Vector6r;

	using ::yade::Matrix3r;
	using ::yade::Matrix6r;
	using ::yade::MatrixXr;
	using ::yade::VectorXr;

	using ::yade::AlignedBox2r;
	using ::yade::AlignedBox3r;
	using ::yade::AngleAxisr;
	using ::yade::Quaternionr;

	using ::yade::Matrix3cr;
	using ::yade::Matrix6cr;
	using ::yade::MatrixXcr;
	using ::yade::Vector2cr;
	using ::yade::Vector3cr;
	using ::yade::Vector6cr;
	using ::yade::VectorXcr;
}

}

#endif

