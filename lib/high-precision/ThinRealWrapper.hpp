/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// This file is mainly to workaround the fact that boost::python has a built-in converter for long double which is losing 3 digits or precision.
//      see: https://www.boost.org/doc/libs/1_71_0/boost/python/converter/builtin_converters.hpp
//
// This file might be useful later in general if we will need some special traits taken care of for Real type.
// IIRC we had previously a problem with Scalar*Vector, like v2=2*v1; and only v2=v1*2; worked. This can be fixed here too:
//                                https://www.boost.org/doc/libs/1_71_0/libs/utility/operators.htm#two_arg
//
// The implementation is based on https://www.boost.org/doc/libs/1_71_0/libs/utility/operators.htm
//                                https://www.boost.org/doc/libs/1_71_0/libs/utility/operators.htm#example
//                                https://www.boost.org/doc/libs/1_71_0/libs/utility/test/operators_test.cpp
// so that all operators =,+,-,*,/,>,<,==,!= etc on Real are predefined using boost.

#ifndef YADE_REAL_STRONG_TYPEDEF_HPP
#define YADE_REAL_STRONG_TYPEDEF_HPP

#include <boost/config.hpp>
#include <boost/move/traits.hpp>
#include <boost/operators.hpp>
#include <boost/type_traits/has_nothrow_assign.hpp>
#include <boost/type_traits/has_nothrow_constructor.hpp>
#include <boost/type_traits/has_nothrow_copy.hpp>
#include <limits>

/*
(A)
	struct StrongReal : boost::partially_ordered1<                                                                                                         \
	                            StrongReal,                                                                                                                \
	                            boost::field_operators1<                                                                                                   \
	                                    StrongReal,                                                                                                        \
	                                    boost::partially_ordered2<StrongReal, UnderlyingStrongReal, boost::field_operators2<StrongReal, UnderlyingStrongReal>>>> {     \
(B)
	struct StrongReal : boost::operators<StrongReal> {                                                                                                     \
(C)
	class StrongReal : boost::ordered_field_operators2<StrongReal, UnderlyingStrongReal>, boost::ordered_field_operators1<StrongReal> {                    \
(D)
	class StrongReal : boost::ordered_field_operators1<StrongReal> {                                                                                       \
(E) // starts to work
	class StrongReal : boost::ordered_field_operators2<StrongReal, UnderlyingStrongReal>,                                                                  \
	                   boost::ordered_field_operators2<StrongReal, double>,                                                                                \
	                   boost::ordered_field_operators2<StrongReal, long unsigned int>,                                                                     \
	                   boost::ordered_field_operators2<StrongReal, long signed int>,                                                                       \
	                   boost::ordered_field_operators2<StrongReal, unsigned int>,                                                                          \
	                   boost::ordered_field_operators2<StrongReal, signed int>,                                                                            \
	                   boost::ordered_field_operators1<StrongReal> {                                                                                       \
(F)
	class StrongReal : boost::field_operators2<StrongReal, UnderlyingStrongReal>,                                                                          \
	                   boost::field_operators2<StrongReal, double>,                                                                                        \
	                   boost::field_operators2<StrongReal, long unsigned int>,                                                                             \
	                   boost::field_operators2<StrongReal, long signed int>,                                                                               \
	                   boost::field_operators2<StrongReal, unsigned int>,                                                                                  \
	                   boost::field_operators2<StrongReal, signed int>,                                                                                    \
	                   boost::field_operators1<StrongReal>,                                                                                                \
	                   boost::partially_ordered2<StrongReal, UnderlyingStrongReal>,                                                                        \
	                   boost::partially_ordered2<StrongReal, double>,                                                                                      \
	                   boost::partially_ordered2<StrongReal, long unsigned int>,                                                                           \
	                   boost::partially_ordered2<StrongReal, long signed int>,                                                                             \
	                   boost::partially_ordered2<StrongReal, unsigned int>,                                                                                \
	                   boost::partially_ordered2<StrongReal, signed int>,                                                                                  \
	                   boost::partially_ordered1<StrongReal> {                                                                                             \

 */

// According to "Ordering note" if we want to have workinf NaN and Inf, we have to use partially_ordered* operators, not ordered_field_operators*
//    https://www.boost.org/doc/libs/1_71_0/libs/utility/operators.htm#ordering
//
// If we need ordering or field operators for other types than those listed below (double, long int, etc…), just add them.
// I skip float for now. See what will happen.

#define YADE_REAL_STRONG_TYPEDEF(UnderlyingStrongReal, StrongReal)                                                                                             \
	class StrongReal : boost::ordered_field_operators2<StrongReal, UnderlyingStrongReal>,                                                                  \
	                   boost::ordered_field_operators2<StrongReal, double>,                                                                                \
	                   boost::ordered_field_operators2<StrongReal, long unsigned int>,                                                                     \
	                   boost::ordered_field_operators2<StrongReal, long signed int>,                                                                       \
	                   boost::ordered_field_operators2<StrongReal, unsigned int>,                                                                          \
	                   boost::ordered_field_operators2<StrongReal, signed int>,                                                                            \
	                   boost::ordered_field_operators1<StrongReal> {                                                                                       \
	private:                                                                                                                                               \
		UnderlyingStrongReal val;                                                                                                                      \
                                                                                                                                                               \
	public:                                                                                                                                                \
		inline StrongReal(const UnderlyingStrongReal& initVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_copy_constructor<UnderlyingStrongReal>::value)     \
		        : val(initVal)                                                                                                                         \
		{                                                                                                                                              \
		}                                                                                                                                              \
		/* This is for floating point literals without suffix, e.g. 0.5 to be automatically promoted to Real  */                                       \
		template <typename OtherType>                                                                                                                  \
		inline StrongReal(const OtherType& initVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_copy_constructor<UnderlyingStrongReal>::value)                \
		        : val(initVal)                                                                                                                         \
		{                                                                                                                                              \
		}                                                                                                                                              \
		inline StrongReal() BOOST_NOEXCEPT_IF(boost::has_nothrow_default_constructor<UnderlyingStrongReal>::value)                                     \
		        : val()                                                                                                                                \
		{                                                                                                                                              \
		}                                                                                                                                              \
		inline StrongReal(const StrongReal& initVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_copy_constructor<UnderlyingStrongReal>::value)               \
		        : val(initVal.val)                                                                                                                     \
		{                                                                                                                                              \
		}                                                                                                                                              \
		/* assignment constructors */                                                                                                                  \
		inline StrongReal& operator=(const StrongReal& rhs) BOOST_NOEXCEPT_IF(boost::has_nothrow_assign<UnderlyingStrongReal>::value)                  \
		{                                                                                                                                              \
			val = rhs.val;                                                                                                                         \
			return *this;                                                                                                                          \
		}                                                                                                                                              \
		inline StrongReal& operator=(const UnderlyingStrongReal& rhs) BOOST_NOEXCEPT_IF(boost::has_nothrow_assign<UnderlyingStrongReal>::value)        \
		{                                                                                                                                              \
			val = rhs;                                                                                                                             \
			return *this;                                                                                                                          \
		}                                                                                                                                              \
		/* This is for floating point literals without suffix, e.g. 0.5 to be automatically promoted to Real  */                                       \
		template <typename OtherType>                                                                                                                  \
		inline StrongReal& operator=(const OtherType& rhs) BOOST_NOEXCEPT_IF(boost::has_nothrow_assign<UnderlyingStrongReal>::value)                   \
		{                                                                                                                                              \
			val = rhs;                                                                                                                             \
			return *this;                                                                                                                          \
		}                                                                                                                                              \
                                                                                                                                                               \
		/* move constructor */                                                                                                                         \
		inline StrongReal(StrongReal&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<UnderlyingStrongReal>::value)                                \
		        : val(std::move(moveVal.val))                                                                                                          \
		{                                                                                                                                              \
		}                                                                                                                                              \
		inline StrongReal& operator=(StrongReal&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<UnderlyingStrongReal>::value)                     \
		{                                                                                                                                              \
			val = std::move(moveVal.val);                                                                                                          \
			return *this;                                                                                                                          \
		}                                                                                                                                              \
                                                                                                                                                               \
		operator const UnderlyingStrongReal&() const { return val; }                                                                                   \
		operator UnderlyingStrongReal&() { return val; }                                                                                               \
                                                                                                                                                               \
		bool                               operator==(const StrongReal& rhs) const { return val == rhs.val; }                                          \
		bool                               operator<(const StrongReal& rhs) const { return val < rhs.val; }                                            \
		bool                               operator>(const StrongReal& rhs) const { return val > rhs.val; }                                            \
		template <typename OtherType> bool operator==(const OtherType& rhs) const { return val == rhs; }                                               \
		template <typename OtherType> bool operator<(const OtherType& rhs) const { return val < rhs; }                                                 \
		template <typename OtherType> bool operator>(const OtherType& rhs) const { return val > rhs; }                                                 \
                                                                                                                                                               \
		inline StrongReal& operator+=(const StrongReal& x)                                                                                             \
		{                                                                                                                                              \
			val += x.val;                                                                                                                          \
			return *this;                                                                                                                          \
		}                                                                                                                                              \
		template <typename OtherType> inline StrongReal& operator+=(const OtherType& x)                                                                \
		{                                                                                                                                              \
			val += x;                                                                                                                              \
			return *this;                                                                                                                          \
		}                                                                                                                                              \
		inline StrongReal& operator-=(const StrongReal& x)                                                                                             \
		{                                                                                                                                              \
			val -= x.val;                                                                                                                          \
			return *this;                                                                                                                          \
		}                                                                                                                                              \
		template <typename OtherType> inline StrongReal& operator-=(const OtherType& x)                                                                \
		{                                                                                                                                              \
			val -= x;                                                                                                                              \
			return *this;                                                                                                                          \
		}                                                                                                                                              \
		inline StrongReal& operator*=(const StrongReal& x)                                                                                             \
		{                                                                                                                                              \
			val *= x.val;                                                                                                                          \
			return *this;                                                                                                                          \
		}                                                                                                                                              \
		template <typename OtherType> inline StrongReal& operator*=(const OtherType& x)                                                                \
		{                                                                                                                                              \
			val *= x;                                                                                                                              \
			return *this;                                                                                                                          \
		}                                                                                                                                              \
		inline StrongReal& operator/=(const StrongReal& x)                                                                                             \
		{                                                                                                                                              \
			val /= x.val;                                                                                                                          \
			return *this;                                                                                                                          \
		}                                                                                                                                              \
		template <typename OtherType> inline StrongReal& operator/=(const OtherType& x)                                                                \
		{                                                                                                                                              \
			val /= x;                                                                                                                              \
			return *this;                                                                                                                          \
		}                                                                                                                                              \
	};

/* maybe add this:
		StrongReal& operator++()                                                                                                                       \
		{                                                                                                                                              \
			++val;                                                                                                                                 \
			return *this;                                                                                                                          \
		};                                                                                                                                             \
		StrongReal& operator--()                                                                                                                       \
		{                                                                                                                                              \
			--val;                                                                                                                                 \
			return *this;                                                                                                                          \
		};                                                                                                                                             \
 */

namespace boost {
namespace multiprecision {
}
}

#define YADE_REAL_NUMERIC_LIMITS(UnderlyingStrongReal, StrongReal)                                                                                             \
	namespace std {                                                                                                                                        \
		template <> struct numeric_limits<StrongReal> {                                                                                                \
			constexpr static auto is_specialized    = std::numeric_limits<UnderlyingStrongReal>::is_specialized;                                   \
			constexpr static auto is_signed         = std::numeric_limits<UnderlyingStrongReal>::is_signed;                                        \
			constexpr static auto is_integer        = std::numeric_limits<UnderlyingStrongReal>::is_integer;                                       \
			constexpr static auto is_exact          = std::numeric_limits<UnderlyingStrongReal>::is_exact;                                         \
			constexpr static auto has_infinity      = std::numeric_limits<UnderlyingStrongReal>::has_infinity;                                     \
			constexpr static auto has_quiet_NaN     = std::numeric_limits<UnderlyingStrongReal>::has_quiet_NaN;                                    \
			constexpr static auto has_signaling_NaN = std::numeric_limits<UnderlyingStrongReal>::has_signaling_NaN;                                \
			constexpr static auto has_denorm        = std::numeric_limits<UnderlyingStrongReal>::has_denorm;                                       \
			constexpr static auto has_denorm_loss   = std::numeric_limits<UnderlyingStrongReal>::has_denorm_loss;                                  \
			constexpr static auto round_style       = std::numeric_limits<UnderlyingStrongReal>::round_style;                                      \
			constexpr static auto is_iec559         = std::numeric_limits<UnderlyingStrongReal>::is_iec559;                                        \
			constexpr static auto is_bounded        = std::numeric_limits<UnderlyingStrongReal>::is_bounded;                                       \
			constexpr static auto is_modulo         = std::numeric_limits<UnderlyingStrongReal>::is_modulo;                                        \
			constexpr static auto digits            = std::numeric_limits<UnderlyingStrongReal>::digits;                                           \
			constexpr static auto digits10          = std::numeric_limits<UnderlyingStrongReal>::digits10;                                         \
			constexpr static auto max_digits10      = std::numeric_limits<UnderlyingStrongReal>::max_digits10;                                     \
			constexpr static auto radix             = std::numeric_limits<UnderlyingStrongReal>::radix;                                            \
			constexpr static auto min_exponent      = std::numeric_limits<UnderlyingStrongReal>::min_exponent;                                     \
			constexpr static auto min_exponent10    = std::numeric_limits<UnderlyingStrongReal>::min_exponent10;                                   \
			constexpr static auto max_exponent      = std::numeric_limits<UnderlyingStrongReal>::max_exponent;                                     \
			constexpr static auto max_exponent10    = std::numeric_limits<UnderlyingStrongReal>::max_exponent10;                                   \
			constexpr static auto traps             = std::numeric_limits<UnderlyingStrongReal>::traps;                                            \
			constexpr static auto tinyness_before   = std::numeric_limits<UnderlyingStrongReal>::tinyness_before;                                  \
			static inline auto    min() { return static_cast<StrongReal>(std::numeric_limits<UnderlyingStrongReal>::min()); }                      \
			static inline auto    lowest() { return static_cast<StrongReal>(std::numeric_limits<UnderlyingStrongReal>::lowest()); }                \
			static inline auto    max() { return static_cast<StrongReal>(std::numeric_limits<UnderlyingStrongReal>::max()); }                      \
			static inline auto    epsilon() { return static_cast<StrongReal>(std::numeric_limits<UnderlyingStrongReal>::epsilon()); }              \
			static inline auto    round_error() { return static_cast<StrongReal>(std::numeric_limits<UnderlyingStrongReal>::round_error()); }      \
			static inline auto    infinity() { return static_cast<StrongReal>(std::numeric_limits<UnderlyingStrongReal>::infinity()); }            \
			static inline auto    quiet_NaN() { return static_cast<StrongReal>(std::numeric_limits<UnderlyingStrongReal>::quiet_NaN()); }          \
			static inline auto    signaling_NaN() { return static_cast<StrongReal>(std::numeric_limits<UnderlyingStrongReal>::signaling_NaN()); }  \
			static inline auto    denorm_min() { return static_cast<StrongReal>(std::numeric_limits<UnderlyingStrongReal>::denorm_min()); }        \
			/* constexpr static auto float_round_style = std::numeric_limits<UnderlyingStrongReal>::float_round_style ; */                         \
			/* constexpr static auto float_denorm_style= std::numeric_limits<UnderlyingStrongReal>::float_denorm_style; */                         \
		};                                                                                                                                             \
	}

#if 0
inline StrongReal abs(const StrongReal& x)
{
	using namespace boost::multiprecision;
	using namespace std;
	return abs(x);
}
static inline StrongReal sqrt(const StrongReal& a)
{
	using namespace boost::multiprecision;
	using namespace std;
	return StrongReal(sqrt(StrongReal(a)));
};
#endif

#endif

