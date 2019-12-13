/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef YADE_REAL_DETECTION_HPP
#define YADE_REAL_DETECTION_HPP

/* This file takes following defines as inputs:
 *
 * #define YADE_REAL_BIT  int(number of bits)
 * #define YADE_REAL_DEC  int(number of decimal places)
 * #define YADE_REAL_MPFR    // defined when Real is expressed using mpfr (requires mpfr)
 * #define YADE_REAL_BBFLOAT // defined when Real is expressed using boost::multiprecision::cpp_bin_float
 *
 * A hardware-accelerated numerical type is used when available, otherwise mpfr or boost::cpp_bin_float.hpp is used for arbitrary precision. The supported types are following:
 *
 *     type                   bits     decimal places         notes
 *     ---------------------------------------------------------------------------
 *     float,float32_t        32       6                      hardware accelerated
 *     double,float64_t       64       15                     hardware accelerated
 *     long double,float80_t  80       18                     hardware accelerated
 *     float128_t             128      33                     depending on processor type it can be hardware accelerated
 *     mpfr                   Nbit     Nbit/(log(2)/log(10))
 *     boost::cpp_bin_float   Nbit     Nbit/(log(2)/log(10))
 *
 *
 * header <boost/cstdfloat.hpp> guarantees standardized floating-point typedefs having specified widths
 * https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/math_toolkit/specified_typedefs.html
 *
 * Depending on platform the fastest one is chosen.
 * https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/math_toolkit/fastest_typdefs.html
 *
 * TODO: Interval types: https://www.boost.org/doc/libs/1_71_0/libs/multiprecision/doc/html/boost_multiprecision/tut/interval/mpfi.html
 *
 * Some hints how to use float128
 *   https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/math_toolkit/float128_hints.html
 *   https://www.boost.org/doc/libs/1_71_0/libs/multiprecision/doc/html/boost_multiprecision/tut/floats/float128.html
 *   e.g.: "Make sure std lib functions are called unqualified so that the correct overload is found via Argument Dependent Lookup (ADL). So write sqrt(variable) and not std::sqrt(variable)."
 *   Eigen author recommends using boost::multiprecision::float128 ← https://forum.kde.org/viewtopic.php?f=74&t=140253
 *
 * Some hints how to use mpfr
 *   https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/index.html
 *   https://www.boost.org/doc/libs/1_71_0/boost/math/tools/precision.hpp
 *   https://www.boost.org/doc/libs/1_71_0/libs/multiprecision/doc/html/index.html
 *   allocate_stack → faster calculations - allocates numbers on stack instead of on heap, this works only with relatively small precisions.
 *   et_on          → faster calculations, slower compilation → FIXME - compile error in Eigen::abs2(…)
 *   et_off         → slower calculations, faster compilation
 */

#include <boost/cstdfloat.hpp> // Must be the first include https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/math_toolkit/rationale.html

#include <boost/config.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/tools/config.hpp>
#include <cmath>
#include <limits>

#include <Eigen/Core>

/*************************************************************************/
/*************************    float 32 bits     **************************/
/*************************************************************************/
#if YADE_REAL_BIT <= 32
using UnderlyingReal = boost::float_fast32_t;
#define YADE_REAL_MATH_NAMESPACE ::std

/*************************************************************************/
/*************************   double 64 bits     **************************/
/*************************************************************************/
#elif YADE_REAL_BIT <= 64
using UnderlyingReal = boost::float_fast64_t;
#define YADE_REAL_MATH_NAMESPACE ::std

/*************************************************************************/
/************************* long double 80 bits  **************************/
/*************************************************************************/
#elif YADE_REAL_BIT <= 80
using UnderlyingReal = boost::float_fast80_t;
#define YADE_REAL_MATH_NAMESPACE ::std
namespace EigenCostReal {
enum { ReadCost = 1, AddCost = 1, MulCost = 1 };
}

/*************************************************************************/
/*************************  float128 128 bits   **************************/
/*************************************************************************/
#elif YADE_REAL_BIT <= 128
#include <quadmath.h>
#include <boost/multiprecision/float128.hpp>
using UnderlyingReal = boost::multiprecision::float128;
#define YADE_REAL_MATH_NAMESPACE ::boost::multiprecision
namespace EigenCostReal {
enum { ReadCost = 1, AddCost = 2, MulCost = 2 };
}

/*************************************************************************/
/*************************         MPFR         **************************/
/*************************************************************************/
#elif defined(YADE_REAL_MPFR)
#include <boost/multiprecision/mpfr.hpp>
template <unsigned int DecimalPlaces>
using UnderlyingRealBackend = boost::multiprecision::mpfr_float_backend<DecimalPlaces, boost::multiprecision::allocate_stack>;
using UnderlyingReal        = boost::multiprecision::number<UnderlyingRealBackend<YADE_REAL_DEC>, boost::multiprecision::et_off>;
#define YADE_REAL_MATH_NAMESPACE ::boost::multiprecision
namespace EigenCostReal {
enum { ReadCost = Eigen::HugeCost, AddCost = Eigen::HugeCost, MulCost = Eigen::HugeCost };
}

/*************************************************************************/
/************************* boost::cpp_bin_float **************************/
/*************************************************************************/
#elif defined(YADE_REAL_BBFLOAT)
#include <boost/multiprecision/cpp_bin_float.hpp>
template <unsigned int DecimalPlaces> using UnderlyingRealBackend = boost::multiprecision::cpp_bin_float<DecimalPlaces>;
using UnderlyingReal = boost::multiprecision::number<UnderlyingRealBackend<YADE_REAL_DEC>, boost::multiprecision::et_off>;
#define YADE_REAL_MATH_NAMESPACE ::boost::multiprecision
namespace EigenCostReal {
enum { ReadCost = Eigen::HugeCost, AddCost = Eigen::HugeCost, MulCost = Eigen::HugeCost };
}

/*************************************************************************/
#elif defined(YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this)
#include <unsupported/Eigen/MPRealSupport>
using UnderlyingReal = ::mpfr::mpreal;
#define YADE_REAL_MATH_NAMESPACE ::mpfr
#else
#error "Real precision is unspecified, there must be a mistake in CMakeLists.txt, the requested #defines should have been provided."
#endif

/*************************************************************************/
/*************************     AND  FINALLY     **************************/
/*************************     declare Real     **************************/
/*************************************************************************/

#if (YADE_REAL_BIT <= 80) and (YADE_REAL_BIT > 64)
// `long double` needs special consideration to workaround boost::python losing 3 digits precision
#include "ThinRealWrapper.hpp"

namespace yade {
using Real = ::yade::ThinRealWrapper<UnderlyingReal>;
}

static_assert(sizeof(yade::Real)==sizeof(UnderlyingReal));

#include "ThinRealWrapperNumericLimits.hpp"
#else

namespace yade {
using Real = UnderlyingReal;
}

#endif

/*************************************************************************/
/*************************    Math functions    **************************/
/*************************************************************************/
#include "MathFunctions.hpp"
#undef YADE_REAL_MATH_NAMESPACE

/*************************************************************************/
/*************************   Eigen  NumTraits   **************************/
/*************************************************************************/
#if (YADE_REAL_BIT > 64) and (not defined(YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this))
using EigenTraitsReal = yade::Real;
#include "EigenNumTraits.hpp"
#endif

/*************************************************************************/
/************************* Vector3 Matrix3 etc  **************************/
/*************************************************************************/
#include "AllMathTypes.hpp"
#undef YADE_REAL_MATH_NAMESPACE

/*************************************************************************/
/*************************  Some sanity checks  **************************/
/*************************************************************************/
#if defined(YADE_REAL_BBFLOAT) and defined(YADE_REAL_MPFR)
#error "Specify either YADE_REAL_MPFR or YADE_REAL_BBFLOAT"
#endif
#if defined(__INTEL_COMPILER) and (YADE_REAL_BIT > 80) and (YADE_REAL_BIT <= 128)
#warning "Intel compiler notes, see: https://www.boost.org/doc/libs/1_71_0/libs/multiprecision/doc/html/boost_multiprecision/tut/floats/float128.html"
#warning "Intel compiler notes: about using flags -Qoption,cpp,--extended_float_type"
#endif

#endif

