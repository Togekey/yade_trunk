/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// Separate header for all Real IO operations.

// Unfortunately std::to_string cannot be relied upon. This is why: https://en.cppreference.com/w/cpp/string/basic_string/to_string
//
// Hence following the boost multiprecision recommendation we will use stringstream as general conversion method:
//   → https://www.boost.org/doc/libs/1_72_0/doc/html/boost_lexical_cast.html
//
// quote: "For more involved conversions, such as where precision or formatting need tighter control than is offered by
// the default behavior of lexical_cast, the conventional std::stringstream approach is recommended. Where the conversions
// are numeric to numeric, boost::numeric_cast may offer more reasonable behavior than lexical_cast."
//

#pragma once

#include <lib/high-precision/Real.hpp>
#include <sstream>
#include <string>

namespace yade {
namespace math {
	// guaranteed maximum precision
	inline std::string toString(const Real& val)
	{
		// FIXME: similar code is in ToFromPythonConverter.hpp, MathSerialization.hpp, MathFunctions.hpp, extract it to single place.
#ifdef YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this
		const auto digs1 = std::numeric_limits<Real>::digits10 + 1;
#else
                static constexpr auto digs1 = std::numeric_limits<Real>::digits10 + 1;
#endif
		std::ostringstream ss;
		ss << std::setprecision(digs1) << val;
		return ss.str();
	};

	inline std::string toString(const Complex& val)
	{
		std::ostringstream ss;
#ifdef YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this
		const auto digs1 = std::numeric_limits<Real>::digits10 + 1;
#else
                static constexpr auto digs1 = std::numeric_limits<Real>::digits10 + 1;
#endif
		ss << std::setprecision(digs1) << val;
		return ss.str();
	};

	// this is just an inline convenience function. It is the same as using std::stringstream.
	// Here stringstream is constructed directly by calling the function with a `const std::string&` argument.
	inline Real fromString(std::stringstream s)
	{
		Real ret;
		s >> ret;
		return ret;
	};

	inline Complex fromString(std::stringstream s)
	{
		Complex ret;
		s >> ret;
		return ret;
	};

}
using yade::math::fromString;
using yade::math::toString;
}

#if (YADE_REAL_BIT > 64)

namespace std {
inline std::string to_string(const ::yade::math::Real& val) { return ::yade::toString(val); };
inline std::string to_string(const ::yade::math::Complex& val) { return ::yade::toString(val); };
}

#endif

