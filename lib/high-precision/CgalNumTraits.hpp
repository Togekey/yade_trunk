/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#if defined(YADE_CGAL) and (not defined(CGAL_NUM_TRAITS_HPP))
#define CGAL_NUM_TRAITS_HPP

#if (YADE_REAL_BIT > 64)
// when `Real` is not plain old `double` we need to define CGAL NumTraits just like it was the case for EigenNumTraits

#include <lib/high-precision/Real.hpp>
#include <CGAL/Interval_nt.h>
#include <CGAL/NT_converter.h>
#include <CGAL/number_type_basic.h>
// cmake already checks for GMP presence when checking for CGAL
#include <CGAL/gmpxx.h>
// newer CGAL uses boost multiprecision
#include <boost/multiprecision/gmp.hpp>

// The traits are those listed in documentation:
//    https://doc.cgal.org/latest/Algebraic_foundations/index.html
//    https://doc.cgal.org/latest/Algebraic_foundations/group__PkgAlgebraicFoundationsRef.html
//
// The functions used by CGAL are on top of https://doc.cgal.org/latest/Algebraic_foundations/group__PkgAlgebraicFoundationsRef.html
// And are implemented below. If some more are necessary they can be added by redirecting to ::yade::math.
// But then a check for them must be also added to:
//    py/high-precision/_math.cpp
//    py/tests/testMath.py
// To make sure they work correctly for all precision types represented by Real.
//
// There is some extra check #ifdef CGAL_CFG_IEEE_754_BUG for long double case in /usr/include/CGAL/long_double.h.
// So #if defined(CGAL_CFG_IEEE_754_BUG) and (Real == long double) then it falls back to CGAL's long double specific functions.

namespace CGAL {

template <> class Is_valid<::yade::Real> : public CGAL::cpp98::unary_function<::yade::Real, bool> {
public:
	bool operator()(const ::yade::Real& x) const
	{
#if defined(CGAL_CFG_IEEE_754_BUG) and (YADE_REAL_BIT == 80)
		return Is_valid<long double>()(static_cast<long double>(x));
#else
		return (x == x);
#endif
	}
};

template <> class Algebraic_structure_traits<::yade::Real> : public Algebraic_structure_traits_base<::yade::Real, Field_with_kth_root_tag> {
public:
	typedef Tag_false Is_exact;
	typedef Tag_true  Is_numerical_sensitive;

	class Sqrt : public CGAL::cpp98::unary_function<Type, Type> {
	public:
		Type operator()(const Type& x) const
		{
			return ::yade::math::sqrt(x);
			// note: that's what would be called for long double case:
			// return Algebraic_structure_traits<long double>::Sqrt()(static_cast<long double>(x));
		}
	};

	class Kth_root : public CGAL::cpp98::binary_function<int, Type, Type> {
	public:
		Type operator()(int k, const Type& x) const
		{
			CGAL_precondition_msg(k > 0, "'k' must be positive for k-th roots");
			return ::yade::math::pow(x, static_cast<::yade::Real>(1.0) / static_cast<::yade::Real>(k));
			// note: that's what would be called for long double case:
			// return Algebraic_structure_traits<long double>::Kth_root()(k, static_cast<long double>(x));
		};
	};
};


template <> class Real_embeddable_traits<::yade::Real> : public INTERN_RET::Real_embeddable_traits_base<::yade::Real, CGAL::Tag_true> {
public:
	class To_interval : public CGAL::cpp98::unary_function<Type, std::pair<double, double>> {
	public:
		std::pair<double, double> operator()(const Type& x) const { return (Interval_nt<>(static_cast<double>(x)) + Interval_nt<>::smallest()).pair(); }
	};

	class Is_finite : public CGAL::cpp98::unary_function<Type, bool> {
	public:
		bool operator()(const Type& x) const
		{
#if defined(CGAL_CFG_IEEE_754_BUG) and (YADE_REAL_BIT == 80)
			// use CGAL workaround for long double if there is an IEEE compiler bug
			return Real_embeddable_traits<long double>::Is_finite()(static_cast<long double>(x));
#else
			return (x == x) and (is_valid(x - x));
#endif
		}
	};
};

// When faster CGAL computations are needed, we might want to use and specialize converter for /usr/include/CGAL/Lazy_exact_nt.h
template <typename GMP1, typename GMP2>
struct NT_converter<::yade::Real, __gmp_expr<GMP1, GMP2>>
        : public CGAL::cpp98::unary_function<::yade::Real, NT_converter<::yade::Real, __gmp_expr<GMP1, GMP2>>> {
	__gmp_expr<GMP1, GMP2> operator()(const ::yade::Real& a) const
	{
		// this is really slow, but fortunately newer CGAL uses boos::multiprecision, below.
		std::string s = std::to_string(a);
		return __gmp_expr<GMP1, GMP2>(s, 10);
	}
};

template <>
struct NT_converter<::yade::Real, boost::multiprecision::mpq_rational>
        : public CGAL::cpp98::unary_function<::yade::Real, NT_converter<::yade::Real, boost::multiprecision::mpq_rational>> {
	boost::multiprecision::mpq_rational operator()(const ::yade::Real& a) const
	{
#ifdef YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this
		// The non-boost MPFR cannot support fast conversion.
		std::stringstream ss;
		ss << a;
		boost::multiprecision::mpq_rational ret {};
		ss >> ret;
		return ret;
#else
		return boost::multiprecision::mpq_rational(static_cast<::yade::math::UnderlyingReal>(a));
#endif
	}
};

} // namespace CGAL

#endif


#include <CGAL/Filtered_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Triangulation_structural_filtering_traits.h>

namespace CGAL {
// The Exact_predicates_inexact_constructions_kernel used `double`. Make corresponding typedef for Real type.
class EReal : public Filtered_kernel_adaptor<
                      Type_equality_wrapper<Simple_cartesian<::yade::Real>::Base<EReal>::Type, EReal>,
#ifdef CGAL_NO_STATIC_FILTERS
                      false>
#else
                      true>
#endif
{
};

typedef EReal Exact_Real_predicates_inexact_constructions_kernel;

template <> struct Triangulation_structural_filtering_traits<EReal> {
	typedef Tag_true Use_structural_filtering_tag;
};

} // namespace CGAL

#endif

