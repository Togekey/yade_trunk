/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef REAL_TO_FROM_PYTHON_CONVERTER_HPP
#define REAL_TO_FROM_PYTHON_CONVERTER_HPP

#include <boost/python.hpp>

/*************************************************************************/
/*************************        Real          **************************/
/*************************************************************************/

template <typename ArbitraryReal> struct ArbitraryReal_to_python {
	static PyObject* convert(const ArbitraryReal& val)
	{
		std::stringstream ss {};
		// the '+1' is to make sure that there are no conversion errors in the last bit.
		constexpr const auto digs1 = std::numeric_limits<ArbitraryReal>::digits10 + 1;
		ss << std::setprecision(digs1) << val;
		::boost::python::object mpmath = ::boost::python::import("mpmath");
		// http://mpmath.org/doc/current/technical.html
		mpmath.attr("mp").attr("dps")  = int(digs1);
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
	}
};

template <typename T> std::string num_to_string(const T& num, int = 0)
{
	constexpr const auto digs1 = std::numeric_limits<T>::digits10 + 1;
	std::stringstream    ss {};
	if (digs1 <= 16) {
		ss << std::setprecision(digs1) << num;
	} else {
		// The only way to make sure that it is copy-pasteable without loss of precision is to put it inside "…"
		ss << "\"" << std::setprecision(digs1) << num << "\"";
	}
	return ss.str();
}


/*************************************************************************/
/*************************       Complex        **************************/
/*************************************************************************/

template <typename ArbitraryComplex> struct ArbitraryComplex_to_python {
	static PyObject* convert(const ArbitraryComplex& val)
	{
		std::stringstream ss_real {};
		std::stringstream ss_imag {};
		// the '+1' is to make sure that there are no conversion errors in the last bit.
		constexpr const auto digs1 = std::numeric_limits<typename ArbitraryComplex::value_type>::digits10 + 1;
		ss_real << std::setprecision(digs1) << val.real();
		ss_imag << std::setprecision(digs1) << val.imag();
		::boost::python::object mpmath = ::boost::python::import("mpmath");
		// http://mpmath.org/doc/current/technical.html
		mpmath.attr("mp").attr("dps")  = int(digs1);
		::boost::python::object result = mpmath.attr("mpc")(ss_real.str(), ss_imag.str());
		return boost::python::incref(result.ptr());
	}
};

// https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/faq/how_can_i_automatically_convert_.html
template <typename ArbitraryComplex> struct ArbitraryComplex_from_python {
	ArbitraryComplex_from_python() { boost::python::converter::registry::push_back(&convertible, &construct, boost::python::type_id<ArbitraryComplex>()); }
	static void* convertible(PyObject* obj_ptr)
	{
		// only python complex or mpmath.mpc(…) objects are supoprted. Strings are not parsed.
		// However a simple workaround is to write mpmath.mpc("1.211213123123123123123123123","-124234234.111")
		PyComplex_AsCComplex(obj_ptr);
		if (PyErr_Occurred() == nullptr)
			return obj_ptr;
		PyErr_Clear();
		return nullptr;
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		std::istringstream ss_real { ::boost::python::call_method<std::string>(
			::boost::python::expect_non_null(PyObject_GetAttrString(obj_ptr, "real")), "__str__") };
		std::istringstream ss_imag { ::boost::python::call_method<std::string>(
			::boost::python::expect_non_null(PyObject_GetAttrString(obj_ptr, "imag")), "__str__") };
		void*              storage = ((boost::python::converter::rvalue_from_python_storage<ArbitraryComplex>*)(data))->storage.bytes;
		new (storage) ArbitraryComplex;
		ArbitraryComplex*                     val = (ArbitraryComplex*)storage;
		typename ArbitraryComplex::value_type re { 0 }, im { 0 };
		ss_real >> re;
		ss_imag >> im;
		*val              = ArbitraryComplex(re, im); // must explicitly call the constructor, static_cast won't work.
		data->convertible = storage;
	}
};

template <typename T> inline std::string num_to_string(const std::complex<T>& num, int = 0)
{
	constexpr const auto digs1 = std::numeric_limits<T>::digits10 + 1;
	std::string          ret;
	if (num.real() != 0 && num.imag() != 0) {
		if (digs1 <= 16) {
			// don't add "+" in the middle if imag is negative and will start with "-"
			return num_to_string(num.real()) + (num.imag() > 0 ? "+" : "") + num_to_string(num.imag()) + "j";
		} else {
			// make sure it is copy-pasteable without loss of precision
			return "mpc(" + num_to_string(num.real()) + "," + num_to_string(num.imag()) + ")";
		}
	}
	// only imaginary is non-zero: skip the real part, and decrease padding to accomoadate the trailing "j"
	if (num.imag() != 0) {
		if (digs1 <= 16) {
			return num_to_string(num.imag()) + "j";
		} else {
			return "mpc(\"0\"," + num_to_string(num.imag()) + ")";
		}
	}
	if (digs1 <= 16) {
		return num_to_string(num.real());
	} else {
		return "mpc(" + num_to_string(num.real()) + ",\"0\")";
	}
}

#ifdef YADE_THIN_REAL_WRAPPER_HPP
template <> inline std::string num_to_string<::yade::Complex>(const ::yade::Complex& num, int)
{
	return num_to_string(static_cast<std::complex<::yade::math::UnderlyingReal>>(num));
}
#endif

#endif

