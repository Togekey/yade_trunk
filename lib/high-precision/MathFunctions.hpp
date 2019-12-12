/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef YADE_THIN_REAL_WRAPPER_MATH_FUNCIONS_HPP
#define YADE_THIN_REAL_WRAPPER_MATH_FUNCIONS_HPP

#include <boost/cstdfloat.hpp> // Must be the first include https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/math_toolkit/rationale.html

#include <boost/config.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/tools/config.hpp>
#include <cmath>
#include <cstdlib>
#include <limits>

// FIXME : use correct perfect forwarding.
//         sometimes I get these errors
//             conversion of first argument to `const Real&&` would be ill-formed
//             cannot bind rvalue reference of type ‘const Real&&’ {aka ‘const ThinRealWrapper<long double>&&’} to lvalue of type ‘ThinRealWrapper<long double>’

#ifdef YADE_THIN_REAL_WRAPPER_HPP
// don't produce errors, even when UnderlyingReal is not related to boost multiprecision
namespace boost {
namespace multiprecision {
}
}

// these are inline redirections towards the correct function
#define YADE_WRAP_FUNC_1(func)                                                                                                                                 \
	inline Real func(const Real& a)                                                                                                                        \
	{                                                                                                                                                      \
		using namespace boost::multiprecision;                                                                                                         \
		using namespace std;                                                                                                                           \
		/*return func(std::forward<UnderlyingReal>(a));*/                                                                                              \
		return ::std::func(static_cast<UnderlyingReal>(a));                                                                                            \
	}

#define YADE_WRAP_FUNC_1_RENAME(func1, func2)                                                                                                                  \
	inline Real func1(const Real& a) { return func2(a); }

#define YADE_WRAP_FUNC_2(func)                                                                                                                                 \
	inline Real func(const Real& a, const Real& b)                                                                                                         \
	{                                                                                                                                                      \
		using namespace boost::multiprecision;                                                                                                         \
		using namespace std;                                                                                                                           \
		return ::std::func(static_cast<UnderlyingReal>(a), static_cast<UnderlyingReal>(b));                                                            \
	}

#define YADE_WRAP_FUNC_2_TYPE1(func, FirstType)                                                                                                                \
	inline Real func(FirstType a, const Real& b)                                                                                                           \
	{                                                                                                                                                      \
		using namespace boost::multiprecision;                                                                                                         \
		using namespace std;                                                                                                                           \
		return ::std::func(a, static_cast<UnderlyingReal>(b));                                                                                         \
	}

#define YADE_WRAP_FUNC_2_TYPE2(func, SecondType)                                                                                                               \
	inline Real func(const Real& a, SecondType b)                                                                                                          \
	{                                                                                                                                                      \
		using namespace boost::multiprecision;                                                                                                         \
		using namespace std;                                                                                                                           \
		return ::std::func(static_cast<UnderlyingReal>(a), b);                                                                                         \
	}

#define YADE_WRAP_FUNC_2_TYPE2_CAST(func, SecondType, CastType)                                                                                                \
	inline Real func(const Real& a, SecondType b)                                                                                                          \
	{                                                                                                                                                      \
		using namespace boost::multiprecision;                                                                                                         \
		using namespace std;                                                                                                                           \
		return ::std::func(static_cast<UnderlyingReal>(a), b->operator CastType());                                                                    \
	}

/*
#define YADE_WRAP_FUNC_2_TYPE2_PTR_TO_PAIR(func, PairType)                                                                                                     \
	inline std::pair<Real, PairType> func(const Real& a)                                                                                                   \
	{                                                                                                                                                      \
		using namespace boost::multiprecision;                                                                                                         \
		using namespace std;                                                                                                                           \
		PairType b;                                                                                                                                    \
		Real     ret = ::std::func(static_cast<UnderlyingReal>(a), &b);                                                                                \
		return std::pair<Real, PairType> { ret, b };                                                                                                   \
	}


#define YADE_WRAP_FUNC_2_TYPE2_PTR_TO_PAIR_CAST(func, PairType, CastType)                                                                                      \
	inline std::pair<Real, PairType> func(const Real& a)                                                                                                   \
	{                                                                                                                                                      \
		using namespace boost::multiprecision;                                                                                                         \
		using namespace std;                                                                                                                           \
		CastType b;                                                                                                                                    \
		Real     ret = ::std::func(static_cast<UnderlyingReal>(a), &b);                                                                                \
		return std::pair<Real, PairType> { ret, b };                                                                                                   \
	}
*/

#define YADE_WRAP_FUNC_3(func)                                                                                                                                 \
	inline Real func(const Real& a, const Real& b, const Real& c)                                                                                          \
	{                                                                                                                                                      \
		using namespace boost::multiprecision;                                                                                                         \
		using namespace std;                                                                                                                           \
		return ::std::func(static_cast<UnderlyingReal>(a), static_cast<UnderlyingReal>(b), static_cast<UnderlyingReal>(c));                            \
	}

#define YADE_WRAP_FUNC_3_TYPE3(func, ThirdType)                                                                                                                \
	inline Real func(const Real& a, const Real& b, ThirdType c)                                                                                            \
	{                                                                                                                                                      \
		using namespace boost::multiprecision;                                                                                                         \
		using namespace std;                                                                                                                           \
		return ::std::func(static_cast<UnderlyingReal>(a), static_cast<UnderlyingReal>(b), c);                                                         \
	}

#else

// TODO: make sure that without ThinRealWrapper the correct functions are called.
#define YADE_WRAP_FUNC_1(func)
#define YADE_WRAP_FUNC_1_RENAME(func1, func2)

#define YADE_WRAP_FUNC_2(func)
#define YADE_WRAP_FUNC_2_TYPE2(func, SecondType, CastType)
#define YADE_WRAP_FUNC_2_TYPE1(func, FirstType)

#define YADE_WRAP_FUNC_3(func)
#define YADE_WRAP_FUNC_3_TYPE3(func, ThirdType)

#endif

// here are listed functions present in
//              https://en.cppreference.com/w/cpp/numeric/math
//              https://en.cppreference.com/w/cpp/numeric/special_functions
// compare with /usr/include/mpreal.h (from debian package libmpfrc++-dev) so we support the same set of functions as MPFR does.

// TODO: https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/math_toolkit/overview_tr1.html
// TODO: They suggest to use this -lboost_math_tr1               boost::math::acosh(x) ↔ boost::math::tr1::acosh(x)
//      ↓ …… for large scale software development where compile times are significant …… difference in performance …… as much as 20 times,
//#include <boost/math/tr1.hpp>

YADE_WRAP_FUNC_1(abs)
YADE_WRAP_FUNC_1(acos)
YADE_WRAP_FUNC_1(acosh)
YADE_WRAP_FUNC_1(asin)
YADE_WRAP_FUNC_1(asinh)
YADE_WRAP_FUNC_1(atan)
YADE_WRAP_FUNC_1(atanh)
YADE_WRAP_FUNC_1(cbrt)
YADE_WRAP_FUNC_1(ceil)
YADE_WRAP_FUNC_1(cos)
YADE_WRAP_FUNC_1(cosh)
YADE_WRAP_FUNC_1(erf)
YADE_WRAP_FUNC_1(erfc)
YADE_WRAP_FUNC_1(exp)
YADE_WRAP_FUNC_1(exp2)
YADE_WRAP_FUNC_1(expm1)
YADE_WRAP_FUNC_1(floor)
YADE_WRAP_FUNC_1(ilogb)
YADE_WRAP_FUNC_1(lgamma)
YADE_WRAP_FUNC_1(log)
YADE_WRAP_FUNC_1(log10)
YADE_WRAP_FUNC_1(log1p)
YADE_WRAP_FUNC_1(log2)
YADE_WRAP_FUNC_1(logb)
//YADE_WRAP_FUNC_1(riemann_zeta) // since C++17
YADE_WRAP_FUNC_1(rint)
YADE_WRAP_FUNC_1(round)
YADE_WRAP_FUNC_1(sin)
YADE_WRAP_FUNC_1(sinh)
YADE_WRAP_FUNC_1(sqrt)
YADE_WRAP_FUNC_1(tan)
YADE_WRAP_FUNC_1(tanh)
YADE_WRAP_FUNC_1(tgamma)
YADE_WRAP_FUNC_1(trunc)

YADE_WRAP_FUNC_1_RENAME(fabs, abs)

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

template <typename T> int sign(T val) { return (T(0) < val) - (val < T(0)); }

YADE_WRAP_FUNC_2(atan2)
//YADE_WRAP_FUNC_2(beta) // since C++17
//YADE_WRAP_FUNC_2(cyl_bessel_i) // since C++17
//YADE_WRAP_FUNC_2(cyl_bessel_j) // since C++17
//YADE_WRAP_FUNC_2(cyl_bessel_k) // since C++17
YADE_WRAP_FUNC_2(fmod)
YADE_WRAP_FUNC_2(hypot)
YADE_WRAP_FUNC_2(max)
YADE_WRAP_FUNC_2(min)
YADE_WRAP_FUNC_2(pow)
YADE_WRAP_FUNC_2(remainder)

//YADE_WRAP_FUNC_2_TYPE1(sph_bessel, unsigned) // since C++17
YADE_WRAP_FUNC_2_TYPE2(ldexp, int)

YADE_WRAP_FUNC_2_TYPE2(frexp, int*) // original C signature
YADE_WRAP_FUNC_2_TYPE2_CAST(modf, Real*, UnderlyingReal*)

//YADE_WRAP_FUNC_2_TYPE2_PTR_TO_PAIR(frexp, int) // for python it returns a pair
//YADE_WRAP_FUNC_2_TYPE2_PTR_TO_PAIR_CAST(modf, Real, UnderlyingReal)

YADE_WRAP_FUNC_3(fma)
//YADE_WRAP_FUNC_3(hypot) // since C++17
YADE_WRAP_FUNC_3_TYPE3(remquo, int*)


#undef YADE_WRAP_FUNC_1
#undef YADE_WRAP_FUNC_1_RENAME

#undef YADE_WRAP_FUNC_2
#undef YADE_WRAP_FUNC_2_TYPE2
#undef YADE_WRAP_FUNC_2_TYPE2_CAST
#undef YADE_WRAP_FUNC_2_TYPE1

#undef YADE_WRAP_FUNC_3
#undef YADE_WRAP_FUNC_3_TYPE3

#endif

