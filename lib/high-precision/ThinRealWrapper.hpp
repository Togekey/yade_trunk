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

	// detect types which are convertible to WrappedReal
	template <typename OtherType> using EnableIfConvertible = std::enable_if_t<std::is_convertible<OtherType, WrappedReal>::value>;

	// detect types which are either WrappedReal or ThinRealWrapper
	// accept all variants: const, &&, const &, etc. In C++20 it will be std::remove_cv_t
	// now it is based on https://en.cppreference.com/w/cpp/types/decay, which may be a little too generous
	// in case of problems we might need to switch to std::remove_cv + std::remove_reference
	template <typename OtherType>
	using EnableIfAnyOfThoseTwo = std::enable_if_t<
	        std::is_same<typename std::decay_t<OtherType>, WrappedReal>::value or std::is_same<typename std::decay_t<OtherType>, ThinRealWrapper>::value>;

	// detect if types are the same, exactly, to the const, volative and &, && qualifiers.
	template <typename OtherType>
	using EnableIfEitherOfThem = std::enable_if_t<std::is_same<OtherType, WrappedReal>::value or std::is_same<OtherType, ThinRealWrapper>::value>;


public:
	// default constructor
	inline ThinRealWrapper() BOOST_NOEXCEPT_IF(boost::has_nothrow_default_constructor<WrappedReal>::value) = default;
	// copy constructor
	inline ThinRealWrapper(const ThinRealWrapper& initVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_copy_constructor<WrappedReal>::value) = default;
	// copy assignment operator
	inline ThinRealWrapper& operator=(const ThinRealWrapper& rhs) BOOST_NOEXCEPT_IF(boost::has_nothrow_assign<WrappedReal>::value) = default;
	// move constructor
	inline ThinRealWrapper(ThinRealWrapper&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value) = default;
	// move assignment operator
	inline ThinRealWrapper& operator=(ThinRealWrapper&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value) = default;
	// destructor
	inline ~ThinRealWrapper() noexcept = default;

	/* Note: both of them are implemened below as move/copy constructors.
	// copy constructor from OtherType which is_convertible o WrappedReal
	//
	template <typename OtherType>
	inline ThinRealWrapper(typename boost::enable_if<std::is_convertible<OtherType, WrappedReal>, const OtherType&>::type initVal)
	        BOOST_NOEXCEPT_IF(boost::has_nothrow_copy_constructor<WrappedReal>::value)
	        : val(initVal)
	{
	}
	// copy assignment operator from OtherType which is_convertible o WrappedReal
	template <typename OtherType>
	inline ThinRealWrapper& operator=(typename boost::enable_if<std::is_convertible<OtherType, WrappedReal>, const OtherType&>::type rhs)
	        BOOST_NOEXCEPT_IF(boost::has_nothrow_assign<WrappedReal>::value)
	{
		val = rhs;
		return *this;
	}
*/
	// move constructor from WrappedReal
	inline ThinRealWrapper(WrappedReal&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value)
	        : val(std::move(moveVal))
	{
	}
	// move constructor from const WrappedReal - is it necessary?
	inline ThinRealWrapper(const WrappedReal&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value)
	        : val(std::move(moveVal))
	{
	}

	// move/copy constructor from OtherType which is_convertible o WrappedReal
	template <typename OtherType, typename = EnableIfConvertible<OtherType>>
	inline ThinRealWrapper(OtherType&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value)
	        : val(std::forward<OtherType>(moveVal))
	{
	}

	// move assignment from WrappedReal
	inline ThinRealWrapper& operator=(WrappedReal&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value)
	{
		val = std::move(moveVal);
		return *this;
	}

	// move assignment from const WrappedReal - is it necessary?
	inline ThinRealWrapper& operator=(const WrappedReal&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value)
	{
		val = std::move(moveVal);
		return *this;
	}

	// move/copy assignment from OtherType which is_convertible o WrappedReal
	template <typename OtherType, typename = EnableIfConvertible<OtherType>>
	inline ThinRealWrapper& operator=(OtherType&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value)
	{
		val = std::forward<OtherType>(moveVal);
		return *this;
	}

	/** FIXME - this should work */
	// perfect forwarding / accessors
	//	template <typename OtherType, typename = EnableIfAnyOfThoseTwo<OtherType>> operator OtherType() { return std::forward<OtherType>(val); }
	//	template <typename OtherType, typename = EnableIfAnyOfThoseTwo<OtherType>> operator const OtherType() const { return std::forward<OtherType>(val); }


	/** FIXME - this almost works */
	// perfect forwarding
	template <typename OtherType, typename = EnableIfEitherOfThem<OtherType>> operator OtherType &&() { return std::forward<OtherType>(val); }
	template <typename OtherType, typename = EnableIfEitherOfThem<OtherType>> operator const OtherType &&() { return std::forward<OtherType>(val); }
	// accessors
	operator const ThinRealWrapper&() const { return val; }
	operator ThinRealWrapper&() { return val; }
	operator const WrappedReal&() const { return val; }
	operator WrappedReal&() { return val; }


	/** FIXME - these are experiments */
	//template <typename OtherType, typename = EnableIfEitherOfThem<OtherType>> operator const OtherType&() const { return val; }
	//template <typename OtherType, typename = EnableIfEitherOfThem<OtherType>> operator OtherType&() { return val; }
	//template <typename OtherType, typename = EnableIfEitherOfThem<OtherType>> operator OtherType() const { return val; }
	//template <typename OtherType, typename = EnableIfEitherOfThem<OtherType>> operator OtherType() { return val; }

	// perfect forwarding
	//	operator ThinRealWrapper &&() { return std::move(*this); }
	//	operator WrappedReal &&() { return std::move(val); }
	//template <typename OtherType, typename = EnableIfAnyOfThoseTwo<OtherType>> operator OtherType &&() { return std::forward<OtherType>(val); }
	//template <typename OtherType, typename = EnableIfAnyOfThoseTwo<OtherType>> operator const OtherType &&() { return std::forward<OtherType>(val); }

	// accessors
	//template <typename OtherType, typename = EnableIfAnyOfThoseTwo<OtherType>> operator OtherType() { return val; }
	//template <typename OtherType, typename = EnableIfAnyOfThoseTwo<OtherType>> operator const OtherType() const { return val; }

	//	operator const ThinRealWrapper&() const { return val; }
	//	operator ThinRealWrapper&() { return val; }
	//	operator const WrappedReal&() const { return val; }
	//	operator WrappedReal&() { return val; }
	// duplicate accessors ?
	//template <typename OtherType, typename = EnableIfAnyOfThoseTwo<OtherType>> operator const OtherType&() const { return val; }
	//template <typename OtherType, typename = EnableIfAnyOfThoseTwo<OtherType>> operator OtherType&() { return val; }

	//template <typename OtherType, typename = EnableIfAnyOfThoseTwo<OtherType>> operator OtherType() const { return val; }
	//template <typename OtherType, typename = EnableIfAnyOfThoseTwo<OtherType>> operator OtherType() { return val; }
};

namespace boost {
namespace multiprecision {
}
}

#define YADE_REAL_NUMERIC_LIMITS(WrappedReal, ThinRealWrapper)                                                                                                   \
	namespace std {                                                                                                                                          \
		template <> struct numeric_limits<ThinRealWrapper> {                                                                                             \
			constexpr static const auto& is_specialized    = std::numeric_limits<WrappedReal>::is_specialized;                                       \
			constexpr static const auto& is_signed         = std::numeric_limits<WrappedReal>::is_signed;                                            \
			constexpr static const auto& is_integer        = std::numeric_limits<WrappedReal>::is_integer;                                           \
			constexpr static const auto& is_exact          = std::numeric_limits<WrappedReal>::is_exact;                                             \
			constexpr static const auto& has_infinity      = std::numeric_limits<WrappedReal>::has_infinity;                                         \
			constexpr static const auto& has_quiet_NaN     = std::numeric_limits<WrappedReal>::has_quiet_NaN;                                        \
			constexpr static const auto& has_signaling_NaN = std::numeric_limits<WrappedReal>::has_signaling_NaN;                                    \
			constexpr static const auto& has_denorm        = std::numeric_limits<WrappedReal>::has_denorm;                                           \
			constexpr static const auto& has_denorm_loss   = std::numeric_limits<WrappedReal>::has_denorm_loss;                                      \
			constexpr static const auto& round_style       = std::numeric_limits<WrappedReal>::round_style;                                          \
			constexpr static const auto& is_iec559         = std::numeric_limits<WrappedReal>::is_iec559;                                            \
			constexpr static const auto& is_bounded        = std::numeric_limits<WrappedReal>::is_bounded;                                           \
			constexpr static const auto& is_modulo         = std::numeric_limits<WrappedReal>::is_modulo;                                            \
			constexpr static const auto& digits            = std::numeric_limits<WrappedReal>::digits;                                               \
			constexpr static const auto& digits10          = std::numeric_limits<WrappedReal>::digits10;                                             \
			constexpr static const auto& max_digits10      = std::numeric_limits<WrappedReal>::max_digits10;                                         \
			constexpr static const auto& radix             = std::numeric_limits<WrappedReal>::radix;                                                \
			constexpr static const auto& min_exponent      = std::numeric_limits<WrappedReal>::min_exponent;                                         \
			constexpr static const auto& min_exponent10    = std::numeric_limits<WrappedReal>::min_exponent10;                                       \
			constexpr static const auto& max_exponent      = std::numeric_limits<WrappedReal>::max_exponent;                                         \
			constexpr static const auto& max_exponent10    = std::numeric_limits<WrappedReal>::max_exponent10;                                       \
			constexpr static const auto& traps             = std::numeric_limits<WrappedReal>::traps;                                                \
			constexpr static const auto& tinyness_before   = std::numeric_limits<WrappedReal>::tinyness_before;                                      \
			static inline auto           min() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::min()); }                     \
			static inline auto           lowest() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::lowest()); }               \
			static inline auto           max() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::max()); }                     \
			static inline auto           epsilon() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::epsilon()); }             \
			static inline auto           round_error() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::round_error()); }     \
			static inline auto           infinity() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::infinity()); }           \
			static inline auto           quiet_NaN() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::quiet_NaN()); }         \
			static inline auto           signaling_NaN() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::signaling_NaN()); } \
			static inline auto           denorm_min() { return static_cast<ThinRealWrapper>(std::numeric_limits<WrappedReal>::denorm_min()); }       \
			/* constexpr static auto float_round_style = std::numeric_limits<WrappedReal>::float_round_style ; */                                    \
			/* constexpr static auto float_denorm_style= std::numeric_limits<WrappedReal>::float_denorm_style; */                                    \
		};                                                                                                                                               \
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

