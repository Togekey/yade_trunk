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

#include <boost/math/constants/constants.hpp>
#include <Eigen/Core>

/*************************************************************************/
/*************************        Real          **************************/
/*************************************************************************/

namespace Eigen {
// NOTE: Don't include this file for float, double, long double. Otherwise you will get errors like:
// error: redefinition of ‘struct Eigen::NumTraits<long double>’
// note: previous definition of ‘struct Eigen::NumTraits<long double>’ in /usr/include/eigen3/Eigen/src/Core/NumTraits.h
template <> struct NumTraits<EigenTraitsReal> : GenericNumTraits<EigenTraitsReal> { // NOTE: Don't include this file for float, double, long double.
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
	typedef EigenTraitsReal Nested;

	static constexpr long get_default_prec = std::numeric_limits<UR>::digits;

	// http://cs.swan.ac.uk/~csoliver/ok-sat-library/internet_html/doc/doc/Mpfr/3.0.0/mpfr.html/Exception-Related-Functions.html
	// https://tspiteri.gitlab.io/gmp-mpfr-sys/mpfr/MPFR-Interface.html
	static inline Real highest(long = get_default_prec) { return std::numeric_limits<UR>::max(); }
	static inline Real lowest(long = get_default_prec) { return std::numeric_limits<UR>::lowest(); }

	// Constants
	static inline Real Pi(long = get_default_prec) { return boost::math::constants::pi<UR>(); }
	static inline Real Euler(long = get_default_prec) { return boost::math::constants::euler<UR>(); }
	static inline Real Log2(long = get_default_prec) { return ::yade::log(EigenTraitsReal(2)); }
	static inline Real Catalan(long = get_default_prec) { return boost::math::constants::catalan<UR>(); }

	static inline Real epsilon(long = get_default_prec) { return std::numeric_limits<UR>::epsilon(); }
	static inline Real epsilon(const Real&) { return std::numeric_limits<UR>::epsilon(); }
	static inline Real smallest_positive() { return std::numeric_limits<UR>::min(); }
	static inline int  digits10(long = get_default_prec) { return std::numeric_limits<UR>::digits10; }
	static inline int  digits10(const Real&) { return std::numeric_limits<UR>::digits10; }
	static inline Real dummy_precision() { return epsilon() * ::yade::pow(Real(10), digits10() / Real(10)); }
};

namespace internal {
	template <> inline EigenTraitsReal random<EigenTraitsReal>() { return ::yade::random(); }
	template <> inline EigenTraitsReal random<EigenTraitsReal>(const EigenTraitsReal& a, const EigenTraitsReal& b)
	{
		return a + (b - a) * ::yade::random01();
	}
	inline bool isMuchSmallerThan(const EigenTraitsReal& a, const EigenTraitsReal& b, const EigenTraitsReal& eps)
	{
		return ::yade::abs(a) <= ::yade::abs(b) * eps;
	}
	inline bool isEqualFuzzy(const EigenTraitsReal& a, const EigenTraitsReal& b, const EigenTraitsReal& eps) { return ::yade::abs(a - b) <= eps; }
	inline bool isApprox(const EigenTraitsReal& a, const EigenTraitsReal& b, const EigenTraitsReal& eps) { return isEqualFuzzy(a, b, eps); }
	inline bool isApproxOrLessThan(const EigenTraitsReal& a, const EigenTraitsReal& b, const EigenTraitsReal& eps)
	{
		return a <= b || isEqualFuzzy(a, b, eps);
	}
	template <> inline long double cast<EigenTraitsReal, long double>(const EigenTraitsReal& x) { return (long double)(x); }
	template <> inline double      cast<EigenTraitsReal, double>(const EigenTraitsReal& x) { return double(x); }
	template <> inline long        cast<EigenTraitsReal, long>(const EigenTraitsReal& x) { return long(x); }
	template <> inline int         cast<EigenTraitsReal, int>(const EigenTraitsReal& x) { return int(x); }
} // end namespace internal

/*************************************************************************/
/*************************       Complex        **************************/
/*************************************************************************/

template <> struct NumTraits<EigenTraitsComplex> : GenericNumTraits<EigenTraitsComplex> { // NOTE: Don't include this file for float, double, long double.
	enum { IsInteger             = 0,
	       IsSigned              = 1,
	       IsComplex             = 1,
	       RequireInitialization = 1,
	       ReadCost              = 2 * ::EigenCostReal::ReadCost,
	       AddCost               = 2 * ::EigenCostReal::AddCost,
	       MulCost               = 4 * ::EigenCostReal::MulCost + 2 * ::EigenCostReal::AddCost };

	typedef EigenTraitsComplex::value_type Real;
	typedef EigenTraitsComplex             NonInteger;
	typedef EigenTraitsComplex             Nested;

	static constexpr long get_default_prec = std::numeric_limits<Real>::digits;

	static inline EigenTraitsComplex highest(long = get_default_prec) { return std::numeric_limits<Real>::max(); }
	static inline EigenTraitsComplex lowest(long = get_default_prec) { return std::numeric_limits<Real>::lowest(); }

	// Constants
	static inline EigenTraitsComplex Pi(long = get_default_prec) { return boost::math::constants::pi<Real>(); }
	static inline EigenTraitsComplex Euler(long = get_default_prec) { return boost::math::constants::euler<Real>(); }
	static inline EigenTraitsComplex Log2(long = get_default_prec) { return ::yade::log(Real(2)); }
	static inline EigenTraitsComplex Catalan(long = get_default_prec) { return boost::math::constants::catalan<Real>(); }

	static inline EigenTraitsComplex epsilon(long = get_default_prec) { return std::numeric_limits<Real>::epsilon(); }
	static inline EigenTraitsComplex epsilon(const EigenTraitsComplex&) { return std::numeric_limits<Real>::epsilon(); }
	static inline EigenTraitsComplex smallest_positive() { return std::numeric_limits<Real>::min(); }
	static inline int                digits10(long = get_default_prec) { return std::numeric_limits<Real>::digits10; }
	static inline int                digits10(const EigenTraitsComplex&) { return std::numeric_limits<Real>::digits10; }
	static inline EigenTraitsComplex dummy_precision() { return epsilon() * ::yade::pow(Real(10), digits10() / Real(10)); }
};

/*
namespace internal {
	template <> inline EigenTraitsComplex random<EigenTraitsComplex>() { return ::yade::random(); }
	template <> inline EigenTraitsComplex random<EigenTraitsComplex>(const EigenTraitsComplex& a, const EigenTraitsComplex& b)
	{
		return a + (b - a) * ::yade::random01();
	}
	inline bool isMuchSmallerThan(const EigenTraitsComplex& a, const EigenTraitsComplex& b, const EigenTraitsComplex& eps)
	{
		return ::yade::abs(a) <= ::yade::abs(b) * eps;
	}
	inline bool isEqualFuzzy(const EigenTraitsComplex& a, const EigenTraitsComplex& b, const EigenTraitsComplex& eps) { return ::yade::abs(a - b) <= eps; }
	inline bool isApprox(const EigenTraitsComplex& a, const EigenTraitsComplex& b, const EigenTraitsComplex& eps) { return isEqualFuzzy(a, b, eps); }
	inline bool isApproxOrLessThan(const EigenTraitsComplex& a, const EigenTraitsComplex& b, const EigenTraitsComplex& eps)
	{
		return a <= b || isEqualFuzzy(a, b, eps);
	}
	template <> inline long double cast<EigenTraitsComplex, long double>(const EigenTraitsComplex& x) { return (long double)(x); }
	template <> inline double      cast<EigenTraitsComplex, double>(const EigenTraitsComplex& x) { return double(x); }
	template <> inline long        cast<EigenTraitsComplex, long>(const EigenTraitsComplex& x) { return long(x); }
	template <> inline int         cast<EigenTraitsComplex, int>(const EigenTraitsComplex& x) { return int(x); }
} // end namespace internal
*/

}

/*
//namespace yade::Math { // maybe add this later

inline const EigenTraitsReal& conj(const EigenTraitsReal& x) { return x; }
inline const EigenTraitsReal& real(const EigenTraitsReal& x) { return x; }
inline EigenTraitsReal        imag(const EigenTraitsReal&) { return 0.; }
inline EigenTraitsReal        abs(const EigenTraitsReal& x)
{
	return ::yade::abs(x);
}
inline EigenTraitsReal abs2(const EigenTraitsReal& x) { return x * x; }

//}
*/

#endif

