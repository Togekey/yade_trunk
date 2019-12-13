/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// This file is mainly to workaround the fact that boost::python has a built-in converter for long double which is losing 3 digits of precision.
//      see: https://www.boost.org/doc/libs/1_71_0/boost/python/converter/builtin_converters.hpp
//
// It will be useful later in general if we will need some special traits taken care of for Real type.

#ifndef YADE_THIN_REAL_WRAPPER_HPP
#define YADE_THIN_REAL_WRAPPER_HPP

#include <boost/config.hpp>
#include <boost/move/traits.hpp>
#include <boost/operators.hpp>
#include <boost/type_traits/has_nothrow_assign.hpp>
#include <boost/type_traits/has_nothrow_constructor.hpp>
#include <boost/type_traits/has_nothrow_copy.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>

// The implementation is based on https://www.boost.org/doc/libs/1_71_0/libs/utility/operators.htm
// so that all operators =,+,-,*,/,>,<,==,!= etc on Real are predefined using boost.
// It is tested with https://www.boost.org/doc/libs/1_72_0/libs/math/doc/html/math_toolkit/real_concepts.html
//  #include <boost/concept/assert.hpp>
//  #include <boost/math/concepts/real_type_concept.hpp>
//  int main() { BOOST_CONCEPT_ASSERT((boost::math::concepts::RealTypeConcept<Real>)); }
//
// see also: /usr/include/boost/math/bindings/e_float.hpp
//           /usr/include/mpreal.h
//           https://www.boost.org/doc/libs/1_72_0/libs/math/test/e_float_concept_check.cpp

// it is possible to #define YADE_IGNORE_IEEE_INFINITY_NAN  ← about that see https://www.boost.org/doc/libs/1_71_0/libs/utility/operators.htm#ordering

namespace yade {

template <typename WrappedReal>
#ifdef YADE_IGNORE_IEEE_INFINITY_NAN
class ThinRealWrapper
        : boost::ordered_field_operators1<ThinRealWrapper<WrappedReal>, boost::ordered_field_operators2<ThinRealWrapper<WrappedReal>, WrappedReal>> {
#else
class ThinRealWrapper : boost::partially_ordered1<
                                ThinRealWrapper<WrappedReal>,
                                boost::partially_ordered2<ThinRealWrapper<WrappedReal>, WrappedReal, boost::field_operators1<ThinRealWrapper<WrappedReal>>>> {
#endif
private:
	WrappedReal val;

	// detect types which are convertible to WrappedReal
	template <typename OtherType> using EnableIfConvertible = std::enable_if_t<std::is_convertible<OtherType, WrappedReal>::value>;

public:
	// default constructor
	ThinRealWrapper() BOOST_NOEXCEPT_IF(boost::has_nothrow_default_constructor<WrappedReal>::value) = default;
	// copy constructor
	ThinRealWrapper(const ThinRealWrapper& initVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_copy_constructor<WrappedReal>::value) = default;
	// copy assignment operator
	ThinRealWrapper& operator=(const ThinRealWrapper& rhs) BOOST_NOEXCEPT_IF(boost::has_nothrow_assign<WrappedReal>::value) = default;
	// move constructor
	ThinRealWrapper(ThinRealWrapper&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value) = default;
	// move assignment operator
	ThinRealWrapper& operator=(ThinRealWrapper&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value) = default;
	// destructor
	~ThinRealWrapper() noexcept = default;

	// NOTE: copy and assignment constructors are implemened below as templated move/copy constructors.

	// move/copy constructor from OtherType which is_convertible o WrappedReal
	template <typename OtherType, typename = EnableIfConvertible<OtherType>>
	ThinRealWrapper(OtherType&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value)
	        : val(std::forward<OtherType>(moveVal))
	{
	}

	// move/copy assignment from OtherType which is_convertible o WrappedReal
	template <typename OtherType, typename = EnableIfConvertible<OtherType>>
	ThinRealWrapper& operator=(OtherType&& moveVal) BOOST_NOEXCEPT_IF(boost::has_nothrow_move<WrappedReal>::value)
	{
		val = std::forward<OtherType>(moveVal);
		return *this;
	}

	// conversion operator to other types
	template <typename OtherType, typename = EnableIfConvertible<OtherType>> explicit operator OtherType() const { return static_cast<OtherType>(val); }

	// https://en.cppreference.com/w/cpp/language/cast_operator
	explicit operator WrappedReal*() { return &val; };
	explicit operator const WrappedReal*() const { return &val; };

	// accessors
	operator const ThinRealWrapper&() const { return val; }
	operator ThinRealWrapper&() { return val; }

	// perfect forwarding
	operator ThinRealWrapper &&() { return std::move(*this); }

	// field operators
	ThinRealWrapper& operator+=(const ThinRealWrapper& x)
	{
		val += x.val;
		return *this;
	}
	ThinRealWrapper& operator-=(const ThinRealWrapper& x)
	{
		val -= x.val;
		return *this;
	}
	ThinRealWrapper& operator*=(const ThinRealWrapper& x)
	{
		val *= x.val;
		return *this;
	}
	ThinRealWrapper& operator/=(const ThinRealWrapper& x)
	{
		val /= x.val;
		return *this;
	}
	const ThinRealWrapper  operator-() const { return -val; }
	const ThinRealWrapper& operator+() const { return *this; }

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

	template <typename OtherType, typename = EnableIfConvertible<OtherType>> friend inline bool operator==(OtherType rhs, const ThinRealWrapper& th)
	{
		th.check(rhs);
		return th.val == rhs;
	}
	template <typename OtherType, typename = EnableIfConvertible<OtherType>> friend inline bool operator!=(OtherType rhs, const ThinRealWrapper& th)
	{
		th.check(rhs);
		return th.val != rhs;
	}
#endif

	// Output/ Input
	friend inline std::ostream& operator<<(std::ostream& os, const ThinRealWrapper& v) { return os << v.val; }
	friend inline std::istream& operator>>(std::istream& is, ThinRealWrapper& v)
	{
		is >> v.val;
		return is;
	}
};

}

#endif

