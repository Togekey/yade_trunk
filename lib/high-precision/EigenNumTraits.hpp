/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef YADE_EIGEN_NUM_TRAITS_HPP
#define YADE_EIGEN_NUM_TRAITS_HPP

// compare this file with /usr/include/eigen3/Eigen/src/Core/NumTraits.h
// http://eigen.tuxfamily.org/dox/TopicCustomizing_CustomScalar.html
// http://eigen.tuxfamily.org/dox-3.2/TopicCustomizingEigen.html

#include <Eigen/Core>
#include <boost/math/constants/constants.hpp>
#include <boost/random.hpp>

namespace Eigen {
// NOTE: Don't include this file for float, double, long double. Otherwise you will get errors like:
// error: redefinition of ‘struct Eigen::NumTraits<long double>’
// note: previous definition of ‘struct Eigen::NumTraits<long double>’ in /usr/include/eigen3/Eigen/src/Core/NumTraits.h
template <> struct NumTraits<EigenTraitsReal> : GenericNumTraits<EigenTraitsReal> { // NOTE: Don't include this file for float, double, long double.
	// /\<ReadCost.*=\|\<MulCost.*=\|\<AddCost.*=
	enum { IsInteger             = 0,
	       IsSigned              = 1,
	       IsComplex             = 0,
	       RequireInitialization = 1,
	       ReadCost              = ::EigenCostReal::ReadCost,
	       AddCost               = ::EigenCostReal::AddCost,
	       MulCost               = ::EigenCostReal::MulCost };

	typedef EigenTraitsReal UR;
	typedef EigenTraitsReal Real;
	typedef EigenTraitsReal NonInteger;

	static constexpr long get_default_prec = std::numeric_limits<UR>::digits;

	// http://cs.swan.ac.uk/~csoliver/ok-sat-library/internet_html/doc/doc/Mpfr/3.0.0/mpfr.html/Exception-Related-Functions.html
	// https://tspiteri.gitlab.io/gmp-mpfr-sys/mpfr/MPFR-Interface.html
	static inline Real highest(long = get_default_prec) { return std::numeric_limits<UR>::max(); }
	static inline Real lowest(long = get_default_prec) { return std::numeric_limits<UR>::lowest(); }

	// Constants
	static inline Real Pi(long = get_default_prec) { return boost::math::constants::pi<UR>(); }
	static inline Real Euler(long = get_default_prec) { return boost::math::constants::euler<UR>(); }
	static inline Real Log2(long = get_default_prec)
	{
		using namespace boost::multiprecision;
		using namespace std;
		return log(EigenTraitsReal(2));
	}
	static inline Real Catalan(long = get_default_prec) { return boost::math::constants::catalan<UR>(); }

	static inline Real epsilon(long = get_default_prec) { return std::numeric_limits<UR>::epsilon(); }
	static inline Real epsilon(const Real&) { return std::numeric_limits<UR>::epsilon(); }

	//#ifdef MPREAL_HAVE_DYNAMIC_STD_NUMERIC_LIMITS
	static inline int digits10(long = get_default_prec) { return std::numeric_limits<UR>::digits10; }
	static inline int digits10(const Real&) { return std::numeric_limits<UR>::digits10; }
	//#endif

	static inline Real dummy_precision()
	{
		using namespace boost::multiprecision;
		return epsilon() * pow(Real(10), digits10() / Real(10));
	}
};

namespace internal {
	namespace supportDetail {
		// random number [0,1)
		static inline EigenTraitsReal random01()
		{
			static boost::random::mt19937 gen;
			return boost::random::generate_canonical<EigenTraitsReal, std::numeric_limits<EigenTraitsReal>::digits>(gen);
		}
	}
	template <> inline EigenTraitsReal random<EigenTraitsReal>() { return supportDetail::random01() * 2 - 1; }

	template <> inline EigenTraitsReal random<EigenTraitsReal>(const EigenTraitsReal& a, const EigenTraitsReal& b)
	{
		return a + (b - a) * supportDetail::random01();
	}

	inline bool isMuchSmallerThan(const EigenTraitsReal& a, const EigenTraitsReal& b, const EigenTraitsReal& eps)
	{
		using namespace boost::multiprecision;
		using namespace std;
		return abs(a) <= abs(b) * eps;
	}

	inline bool isEqualFuzzy(const EigenTraitsReal& a, const EigenTraitsReal& b, const EigenTraitsReal& eps)
	{
		using namespace boost::multiprecision;
		using namespace std;
		return abs(a - b) <= eps;
	}
	inline bool isApprox(const EigenTraitsReal& a, const EigenTraitsReal& b, const EigenTraitsReal& eps) { return isEqualFuzzy(a, b, eps); }

	inline bool isApproxOrLessThan(const EigenTraitsReal& a, const EigenTraitsReal& b, const EigenTraitsReal& eps)
	{
		return a <= b || isEqualFuzzy(a, b, eps);
	}

	template <> inline long double cast<EigenTraitsReal, long double>(const EigenTraitsReal& x) { return (long double)(x); /* x.toLDouble(); */ }

	template <> inline double cast<EigenTraitsReal, double>(const EigenTraitsReal& x) { return double(x); /* x.toDouble(); */ }

	template <> inline long cast<EigenTraitsReal, long>(const EigenTraitsReal& x) { return long(x); /* x.toLong(); */ }

	template <> inline int cast<EigenTraitsReal, int>(const EigenTraitsReal& x) { return int(x); /* int(x.toLong()); */ }
} // end namespace internal
}

/*
//namespace yade::Math { // maybe add this later

inline const EigenTraitsReal& conj(const EigenTraitsReal& x) { return x; }
inline const EigenTraitsReal& real(const EigenTraitsReal& x) { return x; }
inline EigenTraitsReal        imag(const EigenTraitsReal&) { return 0.; }
inline EigenTraitsReal        abs(const EigenTraitsReal& x)
{
	using namespace boost::multiprecision;
	using namespace std;
	return abs(x);
}
inline EigenTraitsReal abs2(const EigenTraitsReal& x) { return x * x; }

//}
*/

#endif

