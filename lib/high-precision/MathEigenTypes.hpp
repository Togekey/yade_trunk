/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef REAL_EXPOSED_TYPES_HPP
#define REAL_EXPOSED_TYPES_HPP

#include <Eigen/Core>
#include <unsupported/Eigen/AlignedVector3>

// integral type for indices, to avoid compiler warnings with int
typedef Eigen::Matrix<int, 1, 1>::Index Index;
typedef Eigen::Matrix<int, 2, 1>        Vector2i;
typedef Eigen::Matrix<int, 3, 1>        Vector3i;
typedef Eigen::Matrix<int, 6, 1>        Vector6i;

/* exposed types */
typedef Eigen::Matrix<Real, 2, 1>                           Vector2r;
typedef Eigen::Matrix<Real, 3, 1>                           Vector3r;
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

#define _COMPLEX_SUPPORT

#ifdef _COMPLEX_SUPPORT
#include <complex>
using std::complex;
typedef Eigen::Matrix<complex<Real>, 2, 1>                           Vector2cr;
typedef Eigen::Matrix<complex<Real>, 3, 1>                           Vector3cr;
typedef Eigen::Matrix<complex<Real>, 6, 1>                           Vector6cr;
typedef Eigen::Matrix<complex<Real>, Eigen::Dynamic, 1>              VectorXcr;
typedef Eigen::Matrix<complex<Real>, 3, 3>                           Matrix3cr;
typedef Eigen::Matrix<complex<Real>, 6, 6>                           Matrix6cr;
typedef Eigen::Matrix<complex<Real>, Eigen::Dynamic, Eigen::Dynamic> MatrixXcr;
#endif

#endif

