/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// NOTE: add more functions as necessary, but remember to add them in py/high-precision/_math.cpp, py/tests/testMath.py and py/tests/testMathHelper.py

// This file contains mathematical functions available in standard library and boost library.
//     https://en.cppreference.com/w/cpp/numeric/math
//     https://en.cppreference.com/w/cpp/numeric/special_functions
// They have to be provided here as inline redirections towards the correct implementation, depending on what precision type yade is being compiled with.
// This is the only way to make sure that ::std, ::boost::math, ::boost::multiprecision are all called correctly.

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
#include <utility>

#ifndef YADE_REAL_MATH_NAMESPACE
#error "This file cannot be included alone, include Real.hpp instead"
#endif


namespace yade {
namespace math {
	/********************************************************************************************/
	/**********************            trigonometric functions             **********************/
	/********************************************************************************************/
	inline Real sin(const Real& a)
	{
		using ::std::sin;
		using YADE_REAL_MATH_NAMESPACE::sin;
		return sin(static_cast<const UnderlyingReal&>(a));
	}
	inline Real sinh(const Real& a)
	{
		using ::std::sinh;
		using YADE_REAL_MATH_NAMESPACE::sinh;
		return sinh(static_cast<const UnderlyingReal&>(a));
	}
	inline Real cos(const Real& a)
	{
		using ::std::cos;
		using YADE_REAL_MATH_NAMESPACE::cos;
		return cos(static_cast<const UnderlyingReal&>(a));
	}
	inline Real cosh(const Real& a)
	{
		using ::std::cosh;
		using YADE_REAL_MATH_NAMESPACE::cosh;
		return cosh(static_cast<const UnderlyingReal&>(a));
	}
	inline Real tan(const Real& a)
	{
		using ::std::tan;
		using YADE_REAL_MATH_NAMESPACE::tan;
		return tan(static_cast<const UnderlyingReal&>(a));
	}
	inline Real tanh(const Real& a)
	{
		using ::std::tanh;
		using YADE_REAL_MATH_NAMESPACE::tanh;
		return tanh(static_cast<const UnderlyingReal&>(a));
	}

	/********************************************************************************************/
	/**********************        complex trigonometric functions         **********************/
	/********************************************************************************************/
	// add more complex functions as necessary, but remember to add them in py/high-precision/_math.cpp, py/tests/testMath.py and py/tests/testMathHelper.py
	inline Complex sin(const Complex& a)
	{
		using ::std::sin;
		using YADE_REAL_MATH_NAMESPACE::sin;
		return sin(static_cast<const std::complex<UnderlyingReal>&>(a));
	}
	inline Complex sinh(const Complex& a)
	{
		using ::std::sinh;
		using YADE_REAL_MATH_NAMESPACE::sinh;
		return sinh(static_cast<const std::complex<UnderlyingReal>&>(a));
	}
	inline Complex cos(const Complex& a)
	{
		using ::std::cos;
		using YADE_REAL_MATH_NAMESPACE::cos;
		return cos(static_cast<const std::complex<UnderlyingReal>&>(a));
	}
	inline Complex cosh(const Complex& a)
	{
		using ::std::cosh;
		using YADE_REAL_MATH_NAMESPACE::cosh;
		return cosh(static_cast<const std::complex<UnderlyingReal>&>(a));
	}
	inline Complex tan(const Complex& a)
	{
		using ::std::tan;
		using YADE_REAL_MATH_NAMESPACE::tan;
		return tan(static_cast<const std::complex<UnderlyingReal>&>(a));
	}
	inline Complex tanh(const Complex& a)
	{
		using ::std::tanh;
		using YADE_REAL_MATH_NAMESPACE::tanh;
		return tanh(static_cast<const std::complex<UnderlyingReal>&>(a));
	}

	/********************************************************************************************/
	/**********************        inverse trigonometric functions         **********************/
	/********************************************************************************************/
	inline Real asin(const Real& a)
	{
		using ::std::asin;
		using YADE_REAL_MATH_NAMESPACE::asin;
		return asin(static_cast<const UnderlyingReal&>(a));
	}
	inline Real asinh(const Real& a)
	{
		using ::std::asinh;
		using YADE_REAL_MATH_NAMESPACE::asinh;
		return 0;//asinh(static_cast<const UnderlyingReal&>(a));
	}
	inline Real acos(const Real& a)
	{
		using ::std::acos;
		using YADE_REAL_MATH_NAMESPACE::acos;
		return acos(static_cast<const UnderlyingReal&>(a));
	}
	inline Real acosh(const Real& a)
	{
		using ::std::acosh;
		using YADE_REAL_MATH_NAMESPACE::acosh;
		return 0;//acosh(static_cast<const UnderlyingReal&>(a));
	}
	inline Real atan(const Real& a)
	{
		using ::std::atan;
		using YADE_REAL_MATH_NAMESPACE::atan;
		return atan(static_cast<const UnderlyingReal&>(a));
	}
	inline Real atanh(const Real& a)
	{
		using ::std::atanh;
		using YADE_REAL_MATH_NAMESPACE::atanh;
		return 0;//atanh(static_cast<const UnderlyingReal&>(a));
	}
	inline Real atan2(const Real& a, const Real& b)
	{
		return YADE_REAL_MATH_NAMESPACE::atan2(static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b));
	}

	/********************************************************************************************/
	/**********************   logarithm, exponential and power functions   **********************/
	/********************************************************************************************/
	inline Real log(const Real& a)
	{
		using ::std::log;
		using YADE_REAL_MATH_NAMESPACE::log;
		return log(static_cast<const UnderlyingReal&>(a));
	}
	inline Real log10(const Real& a)
	{
		using ::std::log10;
		using YADE_REAL_MATH_NAMESPACE::log10;
		return log10(static_cast<const UnderlyingReal&>(a));
	}
	inline Real log1p(const Real& a)
	{
		using ::std::log1p;
		using YADE_REAL_MATH_NAMESPACE::log1p;
		return 0;//log1p(static_cast<const UnderlyingReal&>(a));
	}
	inline Real log2(const Real& a)
	{
		using ::std::log2;
		using YADE_REAL_MATH_NAMESPACE::log2;
		return log2(static_cast<const UnderlyingReal&>(a));
	}
	inline Real logb(const Real& a)
	{
		using ::std::logb;
		using YADE_REAL_MATH_NAMESPACE::logb;
		return logb(static_cast<const UnderlyingReal&>(a));
	}
	inline Real ilogb(const Real& a)
	{
		using ::std::ilogb;
		using YADE_REAL_MATH_NAMESPACE::ilogb;
		return ilogb(static_cast<const UnderlyingReal&>(a));
	}
	inline Real ldexp(const Real& a, int b) { return YADE_REAL_MATH_NAMESPACE::ldexp(static_cast<const UnderlyingReal&>(a), b); }
	// that's original C signature of this function
	inline Real frexp(const Real& a, int* b) { return YADE_REAL_MATH_NAMESPACE::frexp(static_cast<const UnderlyingReal&>(a), b); }
	inline Real exp(const Real& a)
	{
		using ::std::exp;
		using YADE_REAL_MATH_NAMESPACE::exp;
		return exp(static_cast<const UnderlyingReal&>(a));
	}
	inline Real exp2(const Real& a)
	{
		using ::std::exp2;
		using YADE_REAL_MATH_NAMESPACE::exp2;
		return exp2(static_cast<const UnderlyingReal&>(a));
	}
	inline Real expm1(const Real& a)
	{
		using ::std::expm1;
		using YADE_REAL_MATH_NAMESPACE::expm1;
		return 0;//expm1(static_cast<const UnderlyingReal&>(a));
	}
	inline Real pow(const Real& a, const Real& b)
	{
		return YADE_REAL_MATH_NAMESPACE::pow(static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b));
	}
	inline Real sqrt(const Real& a)
	{
		using ::std::sqrt;
		using YADE_REAL_MATH_NAMESPACE::sqrt;
		return sqrt(static_cast<const UnderlyingReal&>(a));
	}
	inline Real cbrt(const Real& a)
	{
		using ::std::cbrt;
		using YADE_REAL_MATH_NAMESPACE::cbrt;
		return 0;//cbrt(static_cast<const UnderlyingReal&>(a));
	}
	inline Real hypot(const Real& a, const Real& b)
	{
		return YADE_REAL_MATH_NAMESPACE::hypot(static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b));
	}
	//YADE_WRAP_FUNC_3(hypot) // since C++17, could be very useful for us

	/********************************************************************************************/
	/**********************        complex logarithm and exponential        *********************/
	/********************************************************************************************/
	// add more complex functions as necessary, but remember to add them in py/high-precision/_math.cpp, py/tests/testMath.py and py/tests/testMathHelper.py
	inline Complex exp(const Complex& a)
	{
		using ::std::exp;
		using YADE_REAL_MATH_NAMESPACE::exp;
		return exp(static_cast<const std::complex<UnderlyingReal>&>(a));
	}
	inline Complex log(const Complex& a)
	{
		using ::std::log;
		using YADE_REAL_MATH_NAMESPACE::log;
		return log(static_cast<const std::complex<UnderlyingReal>&>(a));
	}

	/********************************************************************************************/
	/**********************    min, max, abs, sign, floor, ceil, etc...    **********************/
	/********************************************************************************************/

	// Both must be found by automatic lookup: the ones from ::std and the ones that accept non-double Real types.
	using ::std::abs;
	using ::std::fabs;
	using ::std::max; // this is inside ::yade::math namespace. It is not found by ADL in ::yade namespace when applied to int type or other non-Real type.
	using ::std::min;
#if (defined(YADE_REAL_BIT) and (YADE_REAL_BIT != 64))
	// It turns out that getting min, max to work properly is more tricky than it is for other math functions: https://svn.boost.org/trac10/ticket/11149
	using YADE_REAL_MATH_NAMESPACE::max; // this refers to boost::multiprecision (or eventually to ::mpfr)
	using YADE_REAL_MATH_NAMESPACE::min;
	// make sure that min max can accept (double,Real) argument pairs such as: max(r,0.5);
	inline Real max(const double& a, const Real& b) { return max(static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b)); }
	inline Real min(const double& a, const Real& b) { return min(static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b)); }
	inline Real max(const Real& a, const double& b) { return max(static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b)); }
	inline Real min(const Real& a, const double& b) { return min(static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b)); }
#endif
#if (defined(YADE_REAL_BIT) and (YADE_REAL_BIT > 64))
	inline Real abs(const Real& a)
	{
		using ::std::abs;
		using YADE_REAL_MATH_NAMESPACE::abs;
		return abs(static_cast<const UnderlyingReal&>(a));
	}
	inline Real fabs(const Real& a) { return YADE_REAL_MATH_NAMESPACE::abs(static_cast<const UnderlyingReal&>(a)); }
#endif
	template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }
	template <typename T> int sign(T val) { return (T(0) < val) - (val < T(0)); }

	inline Real floor(const Real& a)
	{
		using ::std::floor;
		using YADE_REAL_MATH_NAMESPACE::floor;
		return floor(static_cast<const UnderlyingReal&>(a));
	}
	inline Real ceil(const Real& a)
	{
		using ::std::ceil;
		using YADE_REAL_MATH_NAMESPACE::ceil;
		return ceil(static_cast<const UnderlyingReal&>(a));
	}
	inline Real round(const Real& a)
	{
		using ::std::round;
		using YADE_REAL_MATH_NAMESPACE::round;
		return round(static_cast<const UnderlyingReal&>(a));
	}
	inline Real rint(const Real& a)
	{
		using ::std::rint;
		using YADE_REAL_MATH_NAMESPACE::rint;
		return rint(static_cast<const UnderlyingReal&>(a));
	}
	inline Real trunc(const Real& a)
	{
		using ::std::trunc;
		using YADE_REAL_MATH_NAMESPACE::trunc;
		return trunc(static_cast<const UnderlyingReal&>(a));
	}

#ifndef YADE_IGNORE_IEEE_INFINITY_NAN
	inline bool isnan(const Real& a)
	{
		using ::std::isnan;
		using YADE_REAL_MATH_NAMESPACE::isnan;
		return isnan(static_cast<const UnderlyingReal&>(a));
	}
	inline bool isinf(const Real& a)
	{
		using ::std::isinf;
		using YADE_REAL_MATH_NAMESPACE::isinf;
		return isinf(static_cast<const UnderlyingReal&>(a));
	}
	inline bool isfinite(const Real& a)
	{
		using ::std::isfinite;
		using YADE_REAL_MATH_NAMESPACE::isfinite;
		return isfinite(static_cast<const UnderlyingReal&>(a));
	}
#endif

	/********************************************************************************************/
	/**********************         complex conj, abs, real, imag          *********************/
	/********************************************************************************************/
	// add more complex functions as necessary, but remember to add them in py/high-precision/_math.cpp and py/tests/testMath.py
	inline Complex conj(const Complex& a) { return ::std::conj(static_cast<const std::complex<UnderlyingReal>&>(a)); }
	inline Real    abs(const Complex& a)
	{
		using ::std::abs;
		using YADE_REAL_MATH_NAMESPACE::abs;
		return abs(static_cast<const std::complex<UnderlyingReal>&>(a));
	}
	inline Real real(const Complex& a) { return ::std::real(static_cast<const std::complex<UnderlyingReal>&>(a)); }
	inline Real imag(const Complex& a) { return ::std::imag(static_cast<const std::complex<UnderlyingReal>&>(a)); }

	/********************************************************************************************/
	/**********************        integer division and remainder          **********************/
	/********************************************************************************************/
	inline Real fmod(const Real& a, const Real& b)
	{
		return YADE_REAL_MATH_NAMESPACE::fmod(static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b));
	}
	inline Real remainder(const Real& a, const Real& b)
	{
		return YADE_REAL_MATH_NAMESPACE::remainder(static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b));
	}
#ifdef YADE_THIN_REAL_WRAPPER_HPP
	inline Real modf(const Real& a, Real* b)
	{
		return YADE_REAL_MATH_NAMESPACE::modf(static_cast<const UnderlyingReal&>(a), b->operator UnderlyingReal*());
	}
#else
	inline Real modf(const Real& a, Real* b) { return YADE_REAL_MATH_NAMESPACE::modf(static_cast<const UnderlyingReal&>(a), b); }
#endif
	inline Real fma(const Real& a, const Real& b, const Real& c)
	{
		return YADE_REAL_MATH_NAMESPACE::fma(
		        static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b), static_cast<const UnderlyingReal&>(c));
	}
	inline Real remquo(const Real& a, const Real& b, int* c)
	{
		return YADE_REAL_MATH_NAMESPACE::remquo(static_cast<const UnderlyingReal&>(a), static_cast<const UnderlyingReal&>(b), c);
	}

	/********************************************************************************************/
	/**********************         special mathematical functions         **********************/
	/********************************************************************************************/
	inline Real erf(const Real& a)
	{
		using ::std::erf;
		using YADE_REAL_MATH_NAMESPACE::erf;
		return 0;//erf(static_cast<const UnderlyingReal&>(a));
	}
	inline Real erfc(const Real& a)
	{
		using ::std::erfc;
		using YADE_REAL_MATH_NAMESPACE::erfc;
		return 0;//erfc(static_cast<const UnderlyingReal&>(a));
	}
	inline Real lgamma(const Real& a)
	{
		using ::std::lgamma;
		using YADE_REAL_MATH_NAMESPACE::lgamma;
		return 0;//lgamma(static_cast<const UnderlyingReal&>(a));
	}

// These will be available in C++17, we could use the ones from boost, if they become necessary.
//YADE_WRAP_FUNC_1(riemann_zeta)
//YADE_WRAP_FUNC_2(beta)
//YADE_WRAP_FUNC_2(cyl_bessel_i)
//YADE_WRAP_FUNC_2(cyl_bessel_j)
//YADE_WRAP_FUNC_2(cyl_bessel_k)
//YADE_WRAP_FUNC_2_TYPE1(sph_bessel, unsigned)


// workaround broken tgamma for boost::float128
#if (defined(YADE_REAL_BIT) and (YADE_REAL_BIT <= 128) and (YADE_REAL_BIT > 80))
	static_assert(std::is_same<UnderlyingReal, boost::multiprecision::float128>::value, "Incorrect type, please file a bug report.");
	inline Real tgamma(const Real& a)
	{
		return 0;
		/*
		if (a >= 0) {
			return YADE_REAL_MATH_NAMESPACE::tgamma(static_cast<UnderlyingReal>(a));
		} else {
			return abs(YADE_REAL_MATH_NAMESPACE::tgamma(static_cast<UnderlyingReal>(a)))
			        * ((static_cast<unsigned long long>(floor(abs(a))) % 2 == 0) ? -1 : 1);
		}
		*/
	}
#else
	inline Real tgamma(const Real& a)
	{
		using ::std::tgamma;
		using YADE_REAL_MATH_NAMESPACE::tgamma;
		return tgamma(static_cast<const UnderlyingReal&>(a));
	}
#endif

	/********************************************************************************************/
	/**********************        extract C-array from std::vector        **********************/
	/********************************************************************************************/

// Some old C library functions need pointer to C-array, this is for compatibility between ThinRealWrapper and UnderlyingReal
#ifdef YADE_THIN_REAL_WRAPPER_HPP
	static_assert(sizeof(Real) == sizeof(UnderlyingReal), "This compiler introduced padding. This breaks binary compatibility");
	static_assert(sizeof(Complex) == sizeof(std::complex<UnderlyingReal>), "This compiler introduced padding, which breaks binary compatibility");

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
		static ::boost::random::mt19937 gen;
		return ::boost::random::generate_canonical<::yade::math::Real, std::numeric_limits<::yade::math::Real>::digits>(gen);
	}

	static inline Real unitRandom() { return random01(); }
	static inline Real random() { return random01() * 2 - 1; }
	static inline Real symmetricRandom() { return random(); }

}
// do we want such alias?
//namespace m = ::yade::math;
//namespace mth = ::yade::math;

}

#endif

