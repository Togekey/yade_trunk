/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef REAL_TO_FROM_PYTHON_CONVERTER_HPP
#define REAL_TO_FROM_PYTHON_CONVERTER_HPP

#include <boost/python.hpp>

#ifdef YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this
#define digits10 digits10()
#endif

#ifdef ARBITRARY_REAL_DEBUG
#include <boost/core/demangle.hpp>
#include <iostream>
template <typename T> std::string infoPrec()
{
	return "\e[93m " + boost::core::demangle(typeid(T).name()) + " prec=" + boost::lexical_cast<std::string>(std::numeric_limits<T>::digits10) + "\e[0m";
}
#endif

template <typename ArbitraryReal> struct ArbitraryReal_to_python {
	static PyObject* convert(const ArbitraryReal& val)
	{
		std::stringstream ss {};
		// the '+1' is to make sure that there are no conversion errors in the last bit.
		ss << std::setprecision(std::numeric_limits<ArbitraryReal>::digits10 + 1) << val;
		::boost::python::object mpmath = ::boost::python::import("mpmath");
#ifdef ARBITRARY_REAL_DEBUG
		std::cerr << "→" << infoPrec<ArbitraryReal>() << "\n"
		          << std::setprecision(std::numeric_limits<ArbitraryReal>::digits10 + 1) << "   HAVE val= " << val << "\n";
		std::cerr << "py::object mpmath pointer is: " << mpmath.ptr() << "\n";
#endif
		// http://mpmath.org/doc/current/technical.html
		mpmath.attr("mp").attr("dps")  = int(std::numeric_limits<ArbitraryReal>::digits10 + 1);
		::boost::python::object result = mpmath.attr("mpf")(ss.str());
		return boost::python::incref(result.ptr());
	}
};

// https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/faq/how_can_i_automatically_convert_.html
template <typename ArbitraryReal> struct ArbitraryReal_from_python {
	ArbitraryReal_from_python() { boost::python::converter::registry::push_back(&convertible, &construct, boost::python::type_id<ArbitraryReal>()); }
	static void* convertible(PyObject* obj_ptr)
	{
		// using long strings or mpmath.mpf(…) object is the only way to get higher precision numbers into C++
		// The line below quickly accepts whatever python is able to convert into float, fortunately this also works for mpmath.mpf(…)
		// this can not work with val=0.123123123123123123123333312312333333123123123, the extra digits are cut-off by python before it reaches this function
		PyFloat_AsDouble(obj_ptr);
		// This quickly returns when argument wasn't a string.
		if (PyErr_Occurred() == nullptr)
			return obj_ptr;
		PyErr_Clear();
		// The quick way didn't work. There was an error, so let's clear it. And check if that is a string with a valid number inside.
		// This is a little more expensive. But it is used very rarely - only when user writes a python line like val="0.123123123123123123123333312312333333123123123"
		// otherwise only mpmath.mpf(NUMBER) objects are passed around inside python scripts which does not reach this line.
		std::istringstream ss { ::boost::python::call_method<std::string>(obj_ptr, "__str__") };
		ArbitraryReal      r;
		ss >> r;
		// Must reach end of string .eof(), otherwise it means there were illegal characters
		return ((not ss.fail()) and (ss.eof())) ? obj_ptr : nullptr;
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		std::istringstream ss { ::boost::python::call_method<std::string>(obj_ptr, "__str__") };

		void* storage = ((boost::python::converter::rvalue_from_python_storage<ArbitraryReal>*)(data))->storage.bytes;
		new (storage) ArbitraryReal;
		ArbitraryReal* val = (ArbitraryReal*)storage;
		ss >> *val;
		data->convertible = storage;

#ifdef ARBITRARY_REAL_DEBUG
		std::cerr << "PyObject* pointer is: " << obj_ptr << " name: " << infoPrec<ArbitraryReal>() << "\n";
		std::cerr << std::setprecision(std::numeric_limits<ArbitraryReal>::digits10 + 1) << "   READ val= " << *val << "\n";
#endif
	}
};

template <typename T> std::string num_to_string(const T& num, int = 0)
{
#ifdef ARBITRARY_REAL_DEBUG
	std::cerr << "\e[91m num_to_string<" << boost::core::demangle(typeid(T).name()) << ">" << (std::numeric_limits<T>::digits10+1) << " number: " << num << "\e[0m\n";
#endif
	std::stringstream ss {};
	ss << std::setprecision(std::numeric_limits<T>::digits10 + 1) << num;
	return ss.str();
}

#ifdef YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this
#undef digits10
#endif

#endif

