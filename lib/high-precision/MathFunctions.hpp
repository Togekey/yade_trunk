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

// NOTE, these using are necessary so that a correct overload will be used regardless of what actually is Real
//		using namespace boost::multiprecision;                                                                                                         \
//		using namespace std;                                                                                                                           \

#define YADE_WRAP_FUNC_1(func)                                                                                                                                 \
	inline const Real func(const Real& a)                                                                                                                  \
	{                                                                                                                                                      \
		using namespace boost::multiprecision;                                                                                                         \
		using namespace std;                                                                                                                           \
		/*return func(std::forward<UnderlyingReal>(a));*/                                                                                              \
		return func(static_cast<UnderlyingReal>(a));                                                                                                   \
		/*return func(a.val);*/                                                                                                                        \
	}

#define YADE_WRAP_FUNC_2(func)                                                                                                                                 \
	inline const Real func(const Real& a, const Real& b)                                                                                                   \
	{                                                                                                                                                      \
		using namespace boost::multiprecision;                                                                                                         \
		using namespace std;                                                                                                                           \
		/*return func(std::forward<UnderlyingReal>(a),std::forward<UnderlyingReal>(b));*/                                                              \
		return func(static_cast<UnderlyingReal>(a), static_cast<UnderlyingReal>(b));                                                                   \
		/*return func(a.val,b.val);*/                                                                                                                  \
	}

// TODO: list all functions present in /usr/include/mpreal.h (from debian package libmpfrc++-dev) so we support the same set of functions as MPFR does.

YADE_WRAP_FUNC_1(log)
YADE_WRAP_FUNC_1(abs)
YADE_WRAP_FUNC_1(sqrt)
YADE_WRAP_FUNC_1(sin)
YADE_WRAP_FUNC_1(cos)
YADE_WRAP_FUNC_1(acos)
YADE_WRAP_FUNC_2(pow)
YADE_WRAP_FUNC_2(atan2)

#endif

