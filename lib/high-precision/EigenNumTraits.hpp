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

/*************************************************************************/
/*************************        Real          **************************/
/*************************************************************************/

namespace Eigen {
// NOTE: Don't include this file for float, double, long double. Otherwise you will get errors like:
// error: redefinition of ‘struct Eigen::NumTraits<long double>’
// note: previous definition of ‘struct Eigen::NumTraits<long double>’ in /usr/include/eigen3/Eigen/src/Core/NumTraits.h
template <> struct NumTraits<::yade::math::Real> : GenericNumTraits<::yade::math::Real> { // NOTE: Don't include this file for float, double, long double.
	enum { IsInteger             = 0,
	       IsSigned              = 1,
	       IsComplex             = 0,
	       RequireInitialization = 1,
	       ReadCost              = ::EigenCostReal::ReadCost,
	       AddCost               = ::EigenCostReal::AddCost,
	       MulCost               = ::EigenCostReal::MulCost };

	typedef ::yade::math::Real Real;
	typedef ::yade::math::Real NonInteger;
	typedef ::yade::math::Real Nested;

	static constexpr long get_default_prec = std::numeric_limits<Real>::digits;

	// http://cs.swan.ac.uk/~csoliver/ok-sat-library/internet_html/doc/doc/Mpfr/3.0.0/mpfr.html/Exception-Related-Functions.html
	// https://tspiteri.gitlab.io/gmp-mpfr-sys/mpfr/MPFR-Interface.html
	static inline Real highest(long = get_default_prec) { return std::numeric_limits<Real>::max(); }
	static inline Real lowest(long = get_default_prec) { return std::numeric_limits<Real>::lowest(); }

	// Constants
	static inline Real Pi(long = get_default_prec) { return boost::math::constants::pi<Real>(); }
	static inline Real Euler(long = get_default_prec) { return boost::math::constants::euler<Real>(); }
	static inline Real Log2(long = get_default_prec) { return ::yade::math::log(Real(2)); }
	static inline Real Catalan(long = get_default_prec) { return boost::math::constants::catalan<Real>(); }

	static inline Real epsilon(long = get_default_prec) { return std::numeric_limits<Real>::epsilon(); }
	static inline Real epsilon(const Real&) { return std::numeric_limits<Real>::epsilon(); }
	static inline Real smallest_positive() { return std::numeric_limits<Real>::min(); }
	static inline int  digits10(long = get_default_prec) { return std::numeric_limits<Real>::digits10; }
	static inline int  digits10(const Real&) { return std::numeric_limits<Real>::digits10; }
	static inline Real dummy_precision() { return epsilon() * ::yade::math::pow(Real(10), digits10() / Real(10)); }
};

namespace internal {
	template <> inline ::yade::math::Real random<::yade::math::Real>() { return ::yade::math::random(); }
	template <> inline ::yade::math::Real random<::yade::math::Real>(const ::yade::math::Real& a, const ::yade::math::Real& b)
	{
		return a + (b - a) * ::yade::math::random01();
	}
	inline bool isMuchSmallerThan(const ::yade::math::Real& a, const ::yade::math::Real& b, const ::yade::math::Real& eps)
	{
		return ::yade::math::abs(a) <= ::yade::math::abs(b) * eps;
	}
	inline bool isEqualFuzzy(const ::yade::math::Real& a, const ::yade::math::Real& b, const ::yade::math::Real& eps)
	{
		return ::yade::math::abs(a - b) <= eps;
	}
	inline bool isApprox(const ::yade::math::Real& a, const ::yade::math::Real& b, const ::yade::math::Real& eps) { return isEqualFuzzy(a, b, eps); }
	inline bool isApproxOrLessThan(const ::yade::math::Real& a, const ::yade::math::Real& b, const ::yade::math::Real& eps)
	{
		return a <= b || isEqualFuzzy(a, b, eps);
	}
	template <> inline long double cast<::yade::math::Real, long double>(const ::yade::math::Real& x) { return (long double)(x); }
	template <> inline double      cast<::yade::math::Real, double>(const ::yade::math::Real& x) { return double(x); }
	template <> inline long        cast<::yade::math::Real, long>(const ::yade::math::Real& x) { return long(x); }
	template <> inline int         cast<::yade::math::Real, int>(const ::yade::math::Real& x) { return int(x); }
} // end namespace internal

/*************************************************************************/
/*************************       Complex        **************************/
/*************************************************************************/

template <> struct NumTraits<::yade::math::Complex> : GenericNumTraits<::yade::math::Complex> { // NOTE: Don't include this file for float, double, long double.
	enum { IsInteger             = 0,
	       IsSigned              = 1,
	       IsComplex             = 1,
	       RequireInitialization = 1,
	       ReadCost              = 2 * ::EigenCostReal::ReadCost,
	       AddCost               = 2 * ::EigenCostReal::AddCost,
	       MulCost               = 4 * ::EigenCostReal::MulCost + 2 * ::EigenCostReal::AddCost };

	typedef ::yade::math::Complex::value_type Real;
	typedef ::yade::math::Complex             Complex;
	typedef ::yade::math::Complex             NonInteger;
	typedef ::yade::math::Complex             Nested;

	static constexpr long get_default_prec = std::numeric_limits<Real>::digits;

	static inline Complex highest(long = get_default_prec) { return std::numeric_limits<Real>::max(); }
	static inline Complex lowest(long = get_default_prec) { return std::numeric_limits<Real>::lowest(); }

	// Constants
	static inline Complex Pi(long = get_default_prec) { return boost::math::constants::pi<Real>(); }
	static inline Complex Euler(long = get_default_prec) { return boost::math::constants::euler<Real>(); }
	static inline Complex Log2(long = get_default_prec) { return ::yade::math::log(Real(2)); }
	static inline Complex Catalan(long = get_default_prec) { return boost::math::constants::catalan<Real>(); }

	static inline Complex epsilon(long = get_default_prec) { return std::numeric_limits<Real>::epsilon(); }
	static inline Complex epsilon(const Complex&) { return std::numeric_limits<Real>::epsilon(); }
	static inline Complex smallest_positive() { return std::numeric_limits<Real>::min(); }
	static inline int     digits10(long = get_default_prec) { return std::numeric_limits<Real>::digits10; }
	static inline int     digits10(const Complex&) { return std::numeric_limits<Real>::digits10; }
	static inline Complex dummy_precision() { return epsilon() * ::yade::math::pow(Real(10), digits10() / Real(10)); }
};

/*
namespace internal {
	template <> inline ::yade::math::Complex random<::yade::math::Complex>() { return ::yade::random(); }
	template <> inline ::yade::math::Complex random<::yade::math::Complex>(const ::yade::math::Complex& a, const ::yade::math::Complex& b)
	{
		return a + (b - a) * ::yade::random01();
	}
	inline bool isMuchSmallerThan(const ::yade::math::Complex& a, const ::yade::math::Complex& b, const ::yade::math::Complex& eps)
	{
		return ::yade::abs(a) <= ::yade::abs(b) * eps;
	}
	inline bool isEqualFuzzy(const ::yade::math::Complex& a, const ::yade::math::Complex& b, const ::yade::math::Complex& eps) { return ::yade::abs(a - b) <= eps; }
	inline bool isApprox(const ::yade::math::Complex& a, const ::yade::math::Complex& b, const ::yade::math::Complex& eps) { return isEqualFuzzy(a, b, eps); }
	inline bool isApproxOrLessThan(const ::yade::math::Complex& a, const ::yade::math::Complex& b, const ::yade::math::Complex& eps)
	{
		return a <= b || isEqualFuzzy(a, b, eps);
	}
	template <> inline long double cast<::yade::math::Complex, long double>(const ::yade::math::Complex& x) { return (long double)(x); }
	template <> inline double      cast<::yade::math::Complex, double>(const ::yade::math::Complex& x) { return double(x); }
	template <> inline long        cast<::yade::math::Complex, long>(const ::yade::math::Complex& x) { return long(x); }
	template <> inline int         cast<::yade::math::Complex, int>(const ::yade::math::Complex& x) { return int(x); }
} // end namespace internal
*/

}

/*
//namespace yade::Math { // maybe add this later

inline const ::yade::math::Real& conj(const ::yade::math::Real& x) { return x; }
inline const ::yade::math::Real& real(const ::yade::math::Real& x) { return x; }
inline ::yade::math::Real        imag(const ::yade::math::Real&) { return 0.; }
inline ::yade::math::Real        abs(const ::yade::math::Real& x)
{
	return ::yade::abs(x);
}
inline ::yade::math::Real abs2(const ::yade::math::Real& x) { return x * x; }

//}
*/

#endif

