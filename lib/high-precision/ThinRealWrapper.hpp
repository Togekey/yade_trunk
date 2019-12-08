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
#include <boost/core/enable_if.hpp>
#include <boost/move/traits.hpp>
#include <boost/operators.hpp>
#include <boost/type_traits/has_nothrow_assign.hpp>
#include <boost/type_traits/has_nothrow_constructor.hpp>
#include <boost/type_traits/has_nothrow_copy.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <limits>
#include <type_traits>

/*
(A)
	struct ThinRealWrapper : boost::partially_ordered1<                                                                                                         \
	                            ThinRealWrapper,                                                                                                                \
	                            boost::field_operators1<                                                                                                   \
	                                    ThinRealWrapper,                                                                                                        \
	                                    boost::partially_ordered2<ThinRealWrapper, WrappedReal, boost::field_operators2<ThinRealWrapper, WrappedReal>>>> {     \
(B)
	struct ThinRealWrapper : boost::operators<ThinRealWrapper> {                                                                                                     \
(C)
	class ThinRealWrapper : boost::ordered_field_operators2<ThinRealWrapper, WrappedReal>, boost::ordered_field_operators1<ThinRealWrapper> {                    \
(D)
	class ThinRealWrapper : boost::ordered_field_operators1<ThinRealWrapper> {                                                                                       \
(E) // starts to work
	class ThinRealWrapper : boost::ordered_field_operators2<ThinRealWrapper, WrappedReal>,                                                                  \
	                   boost::ordered_field_operators2<ThinRealWrapper, double>,                                                                                \
	                   boost::ordered_field_operators2<ThinRealWrapper, long unsigned int>,                                                                     \
	                   boost::ordered_field_operators2<ThinRealWrapper, long signed int>,                                                                       \
	                   boost::ordered_field_operators2<ThinRealWrapper, unsigned int>,                                                                          \
	                   boost::ordered_field_operators2<ThinRealWrapper, signed int>,                                                                            \
	                   boost::ordered_field_operators1<ThinRealWrapper> {                                                                                       \
(F)
	class ThinRealWrapper : boost::field_operators2<ThinRealWrapper, WrappedReal>,                                                                          \
	                   boost::field_operators2<ThinRealWrapper, double>,                                                                                        \
	                   boost::field_operators2<ThinRealWrapper, long unsigned int>,                                                                             \
	                   boost::field_operators2<ThinRealWrapper, long signed int>,                                                                               \
	                   boost::field_operators2<ThinRealWrapper, unsigned int>,                                                                                  \
	                   boost::field_operators2<ThinRealWrapper, signed int>,                                                                                    \
	                   boost::field_operators1<ThinRealWrapper>,                                                                                                \
	                   boost::partially_ordered2<ThinRealWrapper, WrappedReal>,                                                                        \
	                   boost::partially_ordered2<ThinRealWrapper, double>,                                                                                      \
	                   boost::partially_ordered2<ThinRealWrapper, long unsigned int>,                                                                           \
	                   boost::partially_ordered2<ThinRealWrapper, long signed int>,                                                                             \
	                   boost::partially_ordered2<ThinRealWrapper, unsigned int>,                                                                                \
	                   boost::partially_ordered2<ThinRealWrapper, signed int>,                                                                                  \
	                   boost::partially_ordered1<ThinRealWrapper> {                                                                                             \
(G)
	class ThinRealWrapper : boost::ordered_field_operators1<ThinRealWrapper> {                                                                                       \

 */

// According to "Ordering note" if we want to have workinf NaN and Inf, we have to use partially_ordered* operators, not ordered_field_operators*
//    https://www.boost.org/doc/libs/1_71_0/libs/utility/operators.htm#ordering
//
// If we need ordering or field operators for other types than those listed below (double, long int, etc…), just add them.
// I skip float for now. See what will happen.

template <typename WrappedReal> class ThinRealWrapper {
private:
	WrappedReal val;

public:
	// default constructor
	inline ThinRealWrapper() BOOST_NOEXCEPT_IF(boost::has_nothrow_default_constructor<WrappedReal>::value) = default;
	// default constructor
	//inline ThinRealWrapper() BOOST_NOEXCEPT_IF(boost::has_nothrow_default_constructor<WrappedReal>::value)
	//        : val()
	//{
	//}

	// copy constructor
	inline ThinRealWrapper(const ThinRealWrapper& initVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_copy_constructor<WrappedReal>::value) = default;
	// copy constructor
	//inline ThinRealWrapper(const ThinRealWrapper& initVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_copy_constructor<WrappedReal>::value)
	//        : val(initVal.val)
	//{
	//}

	// copy assignment operator
	inline ThinRealWrapper& operator=(const ThinRealWrapper& rhs) BOOST_NOEXCEPT_IF(boost::has_nothrow_assign<WrappedReal>::value) = default;
	// copy assignment operator
	//inline ThinRealWrapper& operator=(const ThinRealWrapper& rhs) BOOST_NOEXCEPT_IF(boost::has_nothrow_assign<WrappedReal>::value)
	//{
	//	val = rhs.val;
	//	return *this;
	//}

	// move constructor
	inline ThinRealWrapper(ThinRealWrapper&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value) = default;
	// move constructor
	//inline ThinRealWrapper(ThinRealWrapper&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value)
	//        : val(std::move(moveVal.val))
	//{
	//}

	// move assignment operator
	inline ThinRealWrapper& operator=(ThinRealWrapper&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value) = default;
	// move assignment operator
	//inline ThinRealWrapper& operator=(ThinRealWrapper&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value)
	//{
	//	val = std::move(moveVal.val);
	//	return *this;
	//}

	// destructor
	inline ~ThinRealWrapper() noexcept = default;

	// copy constructor from OtherType
	template <typename OtherType>
	inline ThinRealWrapper(typename boost::enable_if<std::is_convertible<OtherType, WrappedReal>, const OtherType&>::type initVal)
	        BOOST_NOEXCEPT_IF(boost::has_nothrow_copy_constructor<WrappedReal>::value)
	        : val(initVal)
	{
	}
	// copy assignment operator from OtherType
	template <typename OtherType>
	inline ThinRealWrapper& operator=(typename boost::enable_if<std::is_convertible<OtherType, WrappedReal>, const OtherType&>::type rhs)
	        BOOST_NOEXCEPT_IF(boost::has_nothrow_assign<WrappedReal>::value)
	{
		val = rhs;
		return *this;
	}

	// move constructor from WrappedReal
	inline ThinRealWrapper(WrappedReal&& val) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value)
	        : val(std::move(val))
	{
	}

	inline ThinRealWrapper& operator=(WrappedReal&& val) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value)
	{
		val = std::move(val);
		return *this;
	}

	// accessors
	                              operator const WrappedReal&() const { return val; }
	                              operator WrappedReal&() { return val; }
	template <typename OtherType> operator typename boost::enable_if<std::is_convertible<WrappedReal, OtherType>, OtherType>::type() const { return val; }
	template <typename OtherType> operator typename boost::enable_if<std::is_convertible<WrappedReal, OtherType>, OtherType>::type() { return val; }

	// converters
	//	operator ThinRealWrapper &&() { return std::move(*this); }
	//	operator WrappedReal &&() { return std::move(val); }
};

namespace boost {
namespace multiprecision {
}
}

#define YADE_REAL_NUMERIC_LIMITS(WrappedReal, ThinRealWrapper)                                                                                                 \
	namespace std {                                                                                                                                        \
		template <> struct numeric_limits<ThinRealWrapper> {                                                                                           \
			constexpr static auto is_specialized    = std::numeric_limits<WrappedReal>::is_specialized;                                            \
			constexpr static auto is_signed         = std::numeric_limits<WrappedReal>::is_signed;                                                 \
			constexpr static auto is_integer        = std::numeric_limits<WrappedReal>::is_integer;                                                \
			constexpr static auto is_exact          = std::numeric_limits<WrappedReal>::is_exact;                                                  \
			constexpr static auto has_infinity      = std::numeric_limits<WrappedReal>::has_infinity;                                              \
			constexpr static auto has_quiet_NaN     = std::numeric_limits<WrappedReal>::has_quiet_NaN;                                             \
			constexpr static auto has_signaling_NaN = std::numeric_limits<WrappedReal>::has_signaling_NaN;                                         \
			constexpr static auto has_denorm        = std::numeric_limits<WrappedReal>::has_denorm;                                                \
			constexpr static auto has_denorm_loss   = std::numeric_limits<WrappedReal>::has_denorm_loss;                                           \
			constexpr static auto round_style       = std::numeric_limits<WrappedReal>::round_style;                                               \
			constexpr static auto is_iec559         = std::numeric_limits<WrappedReal>::is_iec559;                                                 \
			constexpr static auto is_bounded        = std::numeric_limits<WrappedReal>::is_bounded;                                                \
			constexpr static auto is_modulo         = std::numeric_limits<WrappedReal>::is_modulo;                                                 \
			constexpr static auto digits            = std::numeric_limits<WrappedReal>::digits;                                                    \
			constexpr static auto digits10          = std::numeric_limits<WrappedReal>::digits10;                                                  \
			constexpr static auto max_digits10      = std::numeric_limits<WrappedReal>::max_digits10;                                              \
			constexpr static auto radix             = std::numeric_limits<WrappedReal>::radix;                                                     \
			constexpr static auto min_exponent      = std::numeric_limits<WrappedReal>::min_exponent;                                              \
			constexpr static auto min_exponent10    = std::numeric_limits<WrappedReal>::min_exponent10;                                            \
			constexpr static auto max_exponent      = std::numeric_limits<WrappedReal>::max_exponent;                                              \
			constexpr static auto max_exponent10    = std::numeric_limits<WrappedReal>::max_exponent10;                                            \
			constexpr static auto traps             = std::numeric_limits<WrappedReal>::traps;                                                     \
			constexpr static auto tinyness_before   = std::numeric_limits<WrappedReal>::tinyness_before;                                           \
			static inline auto    min() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::min()); }                          \
			static inline auto    lowest() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::lowest()); }                    \
			static inline auto    max() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::max()); }                          \
			static inline auto    epsilon() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::epsilon()); }                  \
			static inline auto    round_error() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::round_error()); }          \
			static inline auto    infinity() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::infinity()); }                \
			static inline auto    quiet_NaN() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::quiet_NaN()); }              \
			static inline auto    signaling_NaN() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::signaling_NaN()); }      \
			static inline auto    denorm_min() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::denorm_min()); }            \
			/* constexpr static auto float_round_style = std::numeric_limits<WrappedReal>::float_round_style ; */                                  \
			/* constexpr static auto float_denorm_style= std::numeric_limits<WrappedReal>::float_denorm_style; */                                  \
		};                                                                                                                                             \
	}

#if 0
inline ThinRealWrapper abs(const ThinRealWrapper& x)
{
	using namespace boost::multiprecision;
	using namespace std;
	return abs(x);
}
static inline ThinRealWrapper sqrt(const ThinRealWrapper& a)
{
	using namespace boost::multiprecision;
	using namespace std;
	return ThinRealWrapper(sqrt(ThinRealWrapper(a)));
};
#endif

#endif

