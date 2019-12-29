/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#if defined(YADE_CGAL) and (not defined(CGAL_NUM_TRAITS_HPP))
#define CGAL_NUM_TRAITS_HPP

#include <lib/high-precision/Real.hpp>
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

// inexact predicates do calculations based on inervals.
#if (YADE_REAL_BIT > 64)
std::pair<double, double> to_interval(const ::yade::Real& a) { return ::yade::math::to_interval_target<double>(a); }
double to_double(const ::yade::Real& a) { return a.operator double(); }
#endif

} // namespace CGAL

#endif

