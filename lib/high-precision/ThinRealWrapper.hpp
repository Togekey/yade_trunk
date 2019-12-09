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
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

// Maybe create separate file for ThinRealWrapper IO .hpp ?
#include <iostream>
#include <sstream>

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

(H)
	class ThinRealWrapper : boost::partially_ordered1<ThinRealWrapper<WrappedReal>, boost::field_operators1<ThinRealWrapper<WrappedReal>>> {

 */

// According to "Ordering note" if we want to have workinf NaN and Inf, we have to use partially_ordered* operators, not ordered_field_operators*
//    https://www.boost.org/doc/libs/1_71_0/libs/utility/operators.htm#ordering
//
// If we need ordering or field operators for other types than those listed below (double, long int, etc…), just add them.
// I skip float for now. See what will happen.

//#define YADE_IGNORE_IEEE_INFINITY_NAN

template <typename WrappedReal>
#ifdef YADE_IGNORE_IEEE_INFINITY_NAN
class ThinRealWrapper : boost::ordered_field_operators1<ThinRealWrapper<WrappedReal>> {
#else
class ThinRealWrapper : boost::partially_ordered1<ThinRealWrapper<WrappedReal>, boost::field_operators1<ThinRealWrapper<WrappedReal>>> {
#endif
private:
	WrappedReal val;

	// detect types which are convertible to WrappedReal
	template <typename OtherType> using EnableIfConvertible = std::enable_if_t<std::is_convertible<OtherType, WrappedReal>::value>;


	// detect types which are either WrappedReal or ThinRealWrapper accept all type variants: const, &&, const &, etc. In C++20 it will be std::remove_cv_t
	template <typename OtherType>
	using EnableIfAnyOfThoseTwo = std::enable_if_t<
	        std::is_same<typename std::decay_t<OtherType>, WrappedReal>::value or std::is_same<typename std::decay_t<OtherType>, ThinRealWrapper>::value>;
	// now it is based on https://en.cppreference.com/w/cpp/types/decay, which may be a little too generous
	// in case of problems we might need to switch to std::remove_cv + std::remove_reference


	// detect if types are the same, exactly, to the const, volatile and &, && qualifiers.
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
	// move constructor from const WrappedReal&&  FIXME: is it necessary?
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

	// move assignment from const WrappedReal&&  FIXME: is it necessary?
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


	template <typename OtherType, typename = EnableIfConvertible<OtherType>> explicit operator OtherType() const { return static_cast<OtherType>(val); }
	/*
	// FIXED - replace these with template
	explicit operator bool() const { return val; }
	explicit operator int() const { return val; }
	explicit operator long() const { return val; }
	explicit operator long long() const { return val; }
	explicit operator unsigned() const { return val; }
	explicit operator unsigned long() const { return val; }
	explicit operator unsigned long long() const { return val; }
	explicit operator float() const { return val; }
	explicit operator double() const { return val; }
	explicit operator long double() const { return val; }
*/


	/** FIXME - this should work */
	// accessors / perfect forwarding of &&
	//	template <typename OtherType, typename = EnableIfAnyOfThoseTwo<OtherType>> operator OtherType() { return std::forward<OtherType>(val); }
	//	template <typename OtherType, typename = EnableIfAnyOfThoseTwo<OtherType>> operator const OtherType() const { return std::forward<OtherType>(val); }

	// accessors
	operator const ThinRealWrapper&() const { return val; }
	operator ThinRealWrapper&() { return val; }
	//operator const WrappedReal&() const { return val; }
	//operator WrappedReal&() { return val; }

	// perfect forwarding
	operator ThinRealWrapper &&() { return std::move(*this); }
	//operator WrappedReal &&() { return std::move(val); }
	/** FIXME - this almost works */
	//template <typename OtherType, typename = EnableIfEitherOfThem<OtherType>> operator OtherType &&() { return std::forward<OtherType>(val); }
	//template <typename OtherType, typename = EnableIfEitherOfThem<OtherType>> operator const OtherType &&() { return std::forward<OtherType>(val); }


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

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// ordering operators
	bool operator<(const ThinRealWrapper& rhs) const { return val < rhs.val; }
#ifdef YADE_IGNORE_IEEE_INFINITY_NAN
	bool operator==(const ThinRealWrapper& rhs) const { return val == rhs.val; }
#else
	void check(const ThinRealWrapper& rhs) const
	{
		if (std::isnan(rhs.val) or std::isnan(val) or std::isinf(rhs.val) or std::isinf(val)) {
			throw std::runtime_error("cannot compare NaN, Inf numbers.");
		}
	}
	bool operator==(const ThinRealWrapper& rhs) const
	{
		check(rhs);
		return val == rhs.val;
	}
	bool operator!=(const ThinRealWrapper& rhs) const
	{
		check(rhs);
		return val != rhs.val;
	}
#endif
	//template <typename OtherType> bool operator==(OtherType rhs) const { return val == rhs; }
	//template <typename OtherType> bool operator<(OtherType rhs) const { return val < rhs; }

	//template <typename OtherType, typename = EnableIfConvertible<OtherType>> bool operator==(OtherType&& rhs) const { return val == rhs.val; }
	//template <typename OtherType, typename = EnableIfConvertible<OtherType>> bool operator<(OtherType&& rhs) const { return val < rhs.val; }

	// field operators
	inline ThinRealWrapper& operator+=(const ThinRealWrapper& x)
	{
		val += x.val;
		return *this;
	}
	inline ThinRealWrapper& operator-=(const ThinRealWrapper& x)
	{
		val -= x.val;
		return *this;
	}
	inline ThinRealWrapper& operator*=(const ThinRealWrapper& x)
	{
		val *= x.val;
		return *this;
	}
	inline ThinRealWrapper& operator/=(const ThinRealWrapper& x)
	{
		val /= x.val;
		return *this;
	}
	inline const ThinRealWrapper operator-() const { return -val; }
	/*
	template <typename OtherType> inline ThinRealWrapper& operator+=(const OtherType& x)
	{
		val += x;
		return *this;
	}
	template <typename OtherType> inline ThinRealWrapper& operator-=(const OtherType& x)
	{
		val -= x;
		return *this;
	}
	template <typename OtherType> inline ThinRealWrapper& operator*=(const OtherType& x)
	{
		val *= x;
		return *this;
	}
	template <typename OtherType> inline ThinRealWrapper& operator/=(const OtherType& x)
	{
		val /= x;
		return *this;
	}
*/

	// Output/ Input
	friend inline std::ostream& operator<<(std::ostream& os, const ThinRealWrapper& v) { return os << v.val; }
	friend inline std::istream& operator>>(std::istream& is, ThinRealWrapper& v)
	{
		is >> v.val;
		return is;
	}
};

namespace boost {
namespace multiprecision {
}
}

#endif

