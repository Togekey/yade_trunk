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

template <typename Scalar> using Vector2 = Eigen::Matrix<Scalar, 2, 1>;
template <typename Scalar> using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
template <typename Scalar> using Vector4 = Eigen::Matrix<Scalar, 4, 1>;
template <typename Scalar> using Vector6 = Eigen::Matrix<Scalar, 6, 1>;

template <typename Scalar> using Matrix2 = Eigen::Matrix<Scalar, 3, 3>;
template <typename Scalar> using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
template <typename Scalar> using Matrix4 = Eigen::Matrix<Scalar, 4, 4>;
template <typename Scalar> using Matrix6 = Eigen::Matrix<Scalar, 6, 6>;

/*************************************************************************/
/*************************      Integer         **************************/
/*************************************************************************/
// integral type for indices, to avoid compiler warnings with int
using Index = Eigen::Matrix<int, 1, 1>::Index;

using Vector2i = Vector2<int>;
using Vector3i = Vector3<int>;
using Vector4i = Vector4<int>;
using Vector6i = Vector6<int>;
using VectorXi = Eigen::Matrix<int, Eigen::Dynamic, 1>;

using Matrix2i = Matrix2<int>;
using Matrix3i = Matrix3<int>;
using Matrix4i = Matrix4<int>;
using Matrix6i = Matrix6<int>;
using MatrixXi = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;

/*************************************************************************/
/*************************        Real          **************************/
/*************************************************************************/

// The Vector3r will become AlignedVector3 once we will start using SSE vectorization.
#ifdef EIGEN_DONT_ALIGN
using Vector3r = Vector3<Real>;
#else
using Vector3r = Eigen::AlignedVector3<Real>;
#endif

using Vector2r  = Vector2<Real>;
using Vector3ra = Eigen::AlignedVector3<Real>;
using Vector4r  = Vector4<Real>;
using Vector6r  = Vector6<Real>;
using VectorXr  = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

using Matrix2r = Matrix2<Real>;
using Matrix3r = Matrix3<Real>;
using Matrix4r = Matrix4<Real>;
using Matrix6r = Matrix6<Real>;
using MatrixXr = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

using Quaternionr  = Eigen::Quaternion<Real>;
using AngleAxisr   = Eigen::AngleAxis<Real>;
using AlignedBox3r = Eigen::AlignedBox<Real, 3>;
using AlignedBox2r = Eigen::AlignedBox<Real, 2>;

/*************************************************************************/
/*************************       Complex        **************************/
/*************************************************************************/

using Vector2cr = Vector2<Complex>;
using Vector3cr = Vector3<Complex>;
using Vector4cr = Vector4<Complex>;
using Vector6cr = Vector6<Complex>;
using VectorXcr = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;

using Matrix2cr = Matrix2<Complex>;
using Matrix3cr = Matrix3<Complex>;
using Matrix4cr = Matrix4<Complex>;
using Matrix6cr = Matrix6<Complex>;
using MatrixXcr = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;

/*************************************************************************/
/*************************         Se3          **************************/
/*************************************************************************/

template <class Scalar> class Se3 {
public:
	Vector3<Scalar>           position;
	Eigen::Quaternion<Scalar> orientation;
	Se3() {};
	Se3(Vector3<Scalar> rkP, Eigen::Quaternion<Scalar> qR)
	{
		position    = rkP;
		orientation = qR;
	}
	Se3(const Se3<Scalar>& s)
	{
		position    = s.position;
		orientation = s.orientation;
	}
	Se3(Se3<Scalar>& a, Se3<Scalar>& b)
	{
		position    = b.orientation.inverse() * (a.position - b.position);
		orientation = b.orientation.inverse() * a.orientation;
	}
	Se3<Scalar> inverse() { return Se3(-(orientation.inverse() * position), orientation.inverse()); }
	void        toGLMatrix(float m[16])
	{
		orientation.toGLMatrix(m);
		m[12] = position[0];
		m[13] = position[1];
		m[14] = position[2];
	}
	Vector3<Scalar> operator*(const Vector3<Scalar>& b) { return orientation * b + position; }
	Se3<Scalar>     operator*(const Eigen::Quaternion<Scalar>& b) { return Se3<Scalar>(position, orientation * b); }
	Se3<Scalar>     operator*(const Se3<Scalar>& b) { return Se3<Scalar>(orientation * b.position + position, orientation * b.orientation); }
};

using Se3r = Se3<Real>;

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
	using ::yade::Vector4i;
	using ::yade::Vector4r;
	using ::yade::Vector6i;
	using ::yade::Vector6r;
	using ::yade::VectorXi;
	using ::yade::VectorXr;

	using ::yade::Matrix2i;
	using ::yade::Matrix2r;
	using ::yade::Matrix3i;
	using ::yade::Matrix3r;
	using ::yade::Matrix4i;
	using ::yade::Matrix4r;
	using ::yade::Matrix6i;
	using ::yade::Matrix6r;
	using ::yade::MatrixXi;
	using ::yade::MatrixXr;

	using ::yade::AlignedBox2r;
	using ::yade::AlignedBox3r;
	using ::yade::AngleAxisr;
	using ::yade::Quaternionr;

	using ::yade::Vector2cr;
	using ::yade::Vector3cr;
	using ::yade::Vector4cr;
	using ::yade::Vector6cr;
	using ::yade::VectorXcr;

	using ::yade::Matrix2cr;
	using ::yade::Matrix3cr;
	using ::yade::Matrix4cr;
	using ::yade::Matrix6cr;
	using ::yade::MatrixXcr;

	using ::yade::Se3r;
}

}

#endif

