/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// This file contains mathematical functions available in standard library and boost library.
//     https://en.cppreference.com/w/cpp/numeric/math
//     https://en.cppreference.com/w/cpp/numeric/special_functions
// They have to be provided here as inline redirections towards the correct implementation, depending on what precision type yade is being compiled with.

// TODO: Boost documentation recommends to link with tr1: -lboost_math_tr1 as it provides significant speedup. For example replace boost::math::acosh(x) ↔ boost::math::tr1::acosh(x)
//     https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/math_toolkit/overview_tr1.html
//#include <boost/math/tr1.hpp>


#ifndef YADE_THIN_REAL_WRAPPER_MATH_FUNCIONS_HPP
#define YADE_THIN_REAL_WRAPPER_MATH_FUNCIONS_HPP

#include <boost/config.hpp>
#include <boost/math/complex.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/random.hpp>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <limits>

#ifndef YADE_REAL_MATH_NAMESPACE
#error "This file cannot be included alone, include Real.hpp instead"
#endif

// Macors for quick inline redirections towards the correct function from (1) standard library or (2) boost::multiprecision; depending on which one is used.
#define YADE_WRAP_FUNC_1(func)                                                                                                                                 \
	inline Real func(const Real& a) { return YADE_REAL_MATH_NAMESPACE::func(static_cast<UnderlyingReal>(a)); }

#define YADE_WRAP_FUNC_1_COMPLEX(func)                                                                                                                         \
	inline Complex func(const Complex& a)                                                                                                                  \
	{                                                                                                                                                      \
		using YADE_REAL_MATH_NAMESPACE::func;                                                                                                          \
		using ::std::func;                                                                                                                             \
		return func(static_cast<std::complex<UnderlyingReal>>(a));                                                                                     \
	}

#define YADE_WRAP_FUNC_1_COMPLEX_STD(func)                                                                                                                     \
	inline Complex func(const Complex& a) { return ::std::func(static_cast<std::complex<UnderlyingReal>>(a)); }

#define YADE_WRAP_FUNC_1_COMPLEX_TO_REAL(func)                                                                                                                 \
	inline Real func(const Complex& a)                                                                                                                     \
	{                                                                                                                                                      \
		using YADE_REAL_MATH_NAMESPACE::func;                                                                                                          \
		using ::std::func;                                                                                                                             \
		return func(static_cast<std::complex<UnderlyingReal>>(a));                                                                                     \
	}

#define YADE_WRAP_FUNC_1_COMPLEX_TO_REAL_STD(func)                                                                                                             \
	inline Real func(const Complex& a) { return ::std::func(static_cast<std::complex<UnderlyingReal>>(a)); }

#define YADE_WRAP_FUNC_1_RENAME(func1, func2)                                                                                                                  \
	inline Real func1(const Real& a) { return YADE_REAL_MATH_NAMESPACE::func2(static_cast<UnderlyingReal>(a)); }

#define YADE_WRAP_FUNC_2(func)                                                                                                                                 \
	inline Real func(const Real& a, const Real& b)                                                                                                         \
	{                                                                                                                                                      \
		return YADE_REAL_MATH_NAMESPACE::func(static_cast<UnderlyingReal>(a), static_cast<UnderlyingReal>(b));                                         \
	}

#define YADE_WRAP_FUNC_2_STD(func)                                                                                                                             \
	inline Real func(const Real& a, const Real& b) { return ::std::func(static_cast<UnderlyingReal>(a), static_cast<UnderlyingReal>(b)); }

#define YADE_WRAP_FUNC_2_TYPE1(func, FirstType)                                                                                                                \
	inline Real func(FirstType a, const Real& b) { return YADE_REAL_MATH_NAMESPACE::func(a, static_cast<UnderlyingReal>(b)); }

#define YADE_WRAP_FUNC_2_TYPE2(func, SecondType)                                                                                                               \
	inline Real func(const Real& a, SecondType b) { return YADE_REAL_MATH_NAMESPACE::func(static_cast<UnderlyingReal>(a), b); }

#ifdef YADE_THIN_REAL_WRAPPER_HPP
#define YADE_WRAP_FUNC_2_TYPE2_CAST(func, SecondType, CastType)                                                                                                \
	inline Real func(const Real& a, SecondType b) { return YADE_REAL_MATH_NAMESPACE::func(static_cast<UnderlyingReal>(a), b->operator CastType()); }
#else
#define YADE_WRAP_FUNC_2_TYPE2_CAST(func, SecondType, CastType)                                                                                                \
	inline Real func(const Real& a, SecondType b) { return YADE_REAL_MATH_NAMESPACE::func(static_cast<UnderlyingReal>(a), b); }
#endif

#define YADE_WRAP_FUNC_3(func)                                                                                                                                 \
	inline Real func(const Real& a, const Real& b, const Real& c)                                                                                          \
	{                                                                                                                                                      \
		return YADE_REAL_MATH_NAMESPACE::func(static_cast<UnderlyingReal>(a), static_cast<UnderlyingReal>(b), static_cast<UnderlyingReal>(c));         \
	}

#define YADE_WRAP_FUNC_3_TYPE31(func, ThirdToFirstType)                                                                                                        \
	inline Real func(const Real& a, const Real& b, ThirdToFirstType c)                                                                                     \
	{                                                                                                                                                      \
		return YADE_REAL_MATH_NAMESPACE::func(c, static_cast<UnderlyingReal>(a), static_cast<UnderlyingReal>(b));                                      \
	}

#define YADE_WRAP_FUNC_3_TYPE3(func, ThirdType)                                                                                                                \
	inline Real func(const Real& a, const Real& b, ThirdType c)                                                                                            \
	{                                                                                                                                                      \
		return YADE_REAL_MATH_NAMESPACE::func(static_cast<UnderlyingReal>(a), static_cast<UnderlyingReal>(b), c);                                      \
	}


namespace yade {

/********************************************************************************************/
/**********************            trigonometric functions             **********************/
/********************************************************************************************/

YADE_WRAP_FUNC_1(sin)
YADE_WRAP_FUNC_1(sinh)
YADE_WRAP_FUNC_1(cos)
YADE_WRAP_FUNC_1(cosh)
YADE_WRAP_FUNC_1(tan)
YADE_WRAP_FUNC_1(tanh)

// add more functions as necessary
YADE_WRAP_FUNC_1_COMPLEX(sin)
YADE_WRAP_FUNC_1_COMPLEX(sinh)
YADE_WRAP_FUNC_1_COMPLEX(cos)
YADE_WRAP_FUNC_1_COMPLEX(cosh)
YADE_WRAP_FUNC_1_COMPLEX(tan)
YADE_WRAP_FUNC_1_COMPLEX(tanh)

/********************************************************************************************/
/**********************        inverse trigonometric functions         **********************/
/********************************************************************************************/

YADE_WRAP_FUNC_1(asin)
YADE_WRAP_FUNC_1(asinh)
YADE_WRAP_FUNC_1(acos)
YADE_WRAP_FUNC_1(acosh)
YADE_WRAP_FUNC_1(atan)
YADE_WRAP_FUNC_1(atanh)
YADE_WRAP_FUNC_2(atan2)

/********************************************************************************************/
/**********************   logarithm, exponential and power functions   **********************/
/********************************************************************************************/

YADE_WRAP_FUNC_1(log)
YADE_WRAP_FUNC_1(log10)
YADE_WRAP_FUNC_1(log1p)
YADE_WRAP_FUNC_1(log2)
YADE_WRAP_FUNC_1(logb)
YADE_WRAP_FUNC_1(ilogb)

YADE_WRAP_FUNC_2_TYPE2(ldexp, int)
YADE_WRAP_FUNC_2_TYPE2(frexp, int*) // that's original C signature of this function

YADE_WRAP_FUNC_1(exp)
YADE_WRAP_FUNC_1(exp2)
YADE_WRAP_FUNC_1(expm1)

YADE_WRAP_FUNC_2(pow)
YADE_WRAP_FUNC_1(sqrt)
YADE_WRAP_FUNC_1(cbrt)

YADE_WRAP_FUNC_2(hypot)
//YADE_WRAP_FUNC_3(hypot) // since C++17, could be very useful for us

/********************************************************************************************/
/**********************    min, max, abs, sign, floor, ceil, round     **********************/
/********************************************************************************************/

YADE_WRAP_FUNC_2_STD(max)
YADE_WRAP_FUNC_2_STD(min)
YADE_WRAP_FUNC_1(abs)
YADE_WRAP_FUNC_1_RENAME(fabs, abs)
template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }
template <typename T> int sign(T val) { return (T(0) < val) - (val < T(0)); }
YADE_WRAP_FUNC_1(floor)
YADE_WRAP_FUNC_1(ceil)
YADE_WRAP_FUNC_1(round)
YADE_WRAP_FUNC_1(rint)
YADE_WRAP_FUNC_1(trunc)

YADE_WRAP_FUNC_1_COMPLEX_TO_REAL(abs)

YADE_WRAP_FUNC_1_COMPLEX_STD(conj)
YADE_WRAP_FUNC_1_COMPLEX_TO_REAL_STD(real)
YADE_WRAP_FUNC_1_COMPLEX_TO_REAL_STD(imag)

/********************************************************************************************/
/**********************        integer division and remainder          **********************/
/********************************************************************************************/

YADE_WRAP_FUNC_2(fmod)
YADE_WRAP_FUNC_2(remainder)

#if defined(YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this)
YADE_WRAP_FUNC_2_TYPE2_CAST(modf, Real&, NotUsed)
#else
YADE_WRAP_FUNC_2_TYPE2_CAST(modf, Real*, UnderlyingReal*)
#endif

YADE_WRAP_FUNC_3(fma)
#if defined(YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this)
YADE_WRAP_FUNC_3_TYPE31(remquo, long*)
#else
YADE_WRAP_FUNC_3_TYPE3(remquo, int*)
#endif

/********************************************************************************************/
/**********************         special mathematical functions         **********************/
/********************************************************************************************/

YADE_WRAP_FUNC_1(erf)
YADE_WRAP_FUNC_1(erfc)
YADE_WRAP_FUNC_1(lgamma)
// workaround broken tgamma for boost::float128
#if (YADE_REAL_BIT <= 128) and (YADE_REAL_BIT > 80)
static_assert(std::is_same<UnderlyingReal, boost::multiprecision::float128>::value, "Incorrect type, please file a bug report.");
inline Real tgamma(const Real& a)
{
	if (a >= 0) {
		return YADE_REAL_MATH_NAMESPACE::tgamma(static_cast<UnderlyingReal>(a));
	} else {
		return abs(YADE_REAL_MATH_NAMESPACE::tgamma(static_cast<UnderlyingReal>(a)))
		        * ((static_cast<unsigned long long>(floor(abs(a))) % 2 == 0) ? -1 : 1);
	}
}
#else
YADE_WRAP_FUNC_1(tgamma)
#endif

// These will be available in C++17, we could use the ones from boost, if they become necessary.
//YADE_WRAP_FUNC_1(riemann_zeta)
//YADE_WRAP_FUNC_2(beta)
//YADE_WRAP_FUNC_2(cyl_bessel_i)
//YADE_WRAP_FUNC_2(cyl_bessel_j)
//YADE_WRAP_FUNC_2(cyl_bessel_k)
//YADE_WRAP_FUNC_2_TYPE1(sph_bessel, unsigned)


/********************************************************************************************/
/**********************        extract C-array from std::vector        **********************/
/********************************************************************************************/

// Some old C library functions need pointer to C-array, this is for compatibility between ThinRealWrapper and UnderlyingReal
#ifdef YADE_THIN_REAL_WRAPPER_HPP
static_assert(sizeof(Real) == sizeof(UnderlyingReal), "This compiler introduced padding. This breaks binary compatibility");
static_assert(sizeof(yade::Complex) == sizeof(std::complex<UnderlyingReal>), "This compiler introduced padding, which breaks binary compatibility");

static inline const UnderlyingReal* constVectorData(const std::vector<Real>& v) { return v.data()->operator const UnderlyingReal*(); }
static inline UnderlyingReal*       vectorData(std::vector<Real>& v) { return v.data()->operator UnderlyingReal*(); }
#else
static inline const UnderlyingReal* constVectorData(const std::vector<Real>& v) { return v.data(); }
static inline UnderlyingReal*       vectorData(std::vector<Real>& v) { return v.data(); }
#endif

/********************************************************************************************/
/**********************                     random                     **********************/
/********************************************************************************************/

// These random functions are necessary for Eigen library to for example write in python: Vector3.Random()
// generate random number [0,1)
static inline Real random01()
{
#if defined(YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this)
	return ::mpfr::random();
#else
        static ::boost::random::mt19937 gen;
        return ::boost::random::generate_canonical<Real, std::numeric_limits<Real>::digits>(gen);
#endif
}

static inline Real unitRandom() { return random01(); }
static inline Real random() { return random01() * 2 - 1; }

}

#undef YADE_WRAP_FUNC_1
#undef YADE_WRAP_FUNC_1_RENAME
#undef YADE_WRAP_FUNC_2
#undef YADE_WRAP_FUNC_2_STD
#undef YADE_WRAP_FUNC_2_TYPE2
#undef YADE_WRAP_FUNC_2_TYPE2_CAST
#undef YADE_WRAP_FUNC_2_TYPE1
#undef YADE_WRAP_FUNC_3
#undef YADE_WRAP_FUNC_3_TYPE3
#undef YADE_WRAP_FUNC_3_TYPE31

#endif

