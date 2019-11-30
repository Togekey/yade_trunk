/*************************************************************************
*  2019      Janek Kozicki                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include <core/Omega.hpp>
#include <lib/base/Logging.hpp>
#include <lib/pyutil/doc_opts.hpp>
#include <lib/base/AliasNamespaces.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <limits>
#include <sstream>

#if YADE_REAL_BIT>=128
#include <boost/cstdfloat.hpp>
#include <boost/math/cstdfloat/cstdfloat_types.hpp>
#include <boost/multiprecision/float128.hpp>
#endif

#ifdef YADE_MPFR
#include <boost/multiprecision/mpfr.hpp>
#endif

namespace yade { // Cannot have #include directive inside.

CREATE_CPP_LOCAL_LOGGER("_testCppPy.cpp");

template<typename T>
py::tuple testPi() {
	T ret = boost::math::constants::pi<T>();
	std::stringstream ss;
	ss << std::setprecision(std::numeric_limits<T>::digits10+1) << ret;
	std::cout << __PRETTY_FUNCTION__ << "    \t returns pi = " << ss.str() << "\n";
	return py::make_tuple( ss.str() , ret );
}

} // namespace yade

BOOST_PYTHON_MODULE(_testCppPy) {
	namespace py = ::boost::python;
	YADE_SET_DOCSTRING_OPTS;
// The docstrings syntax https://thomas-cokelaer.info/tutorials/sphinx/docstring_python.html

	py::def("testPiFloat", yade::testPi<float>, R"""(
This function returns a string and ``float`` representation of pi. They can be compared afterwards to see if precision conversion is handled correctly.

:return: ``(string(pi) , float(pi))`` from C++
	)""");

	py::def("testPiDouble", yade::testPi<double>, R"""(
This function returns a string and ``double`` representation of pi. They can be compared afterwards to see if precision conversion is handled correctly.

:return: ``(string(pi) , double(pi))`` from C++
	)""");

	py::def("testPiLongDouble", yade::testPi<long double>, R"""(
This function returns a string and ``long double`` representation of pi. They can be compared afterwards to see if precision conversion is handled correctly.

:return: ``(string(pi) , (long double)(pi))`` from C++
	)""");

#if YADE_REAL_BIT>=128
	py::def("testPiFloat128", yade::testPi<boost::float128_t>, R"""(
This function returns a string and ``boost::float128_t`` representation of pi. They can be compared afterwards to see if precision conversion is handled correctly.

:return: ``(string(pi) , boost::float128_t(pi))`` from C++
	)""");
#endif

#ifdef YADE_MPFR
	#define stringify(ss) #ss   // FIXME : this macro is not working,  FIXME: https://gitlab.com/yade-dev/trunk/issues/97
	py::def("testPiMpfr", yade::testPi<yade::bmp::number<yade::bmp::mpfr_float_backend<YADE_REAL_DEC>>>, R"""(
This function returns a string and ``boost::mpfr_float_backend<)""" stringify(YADE_REAL_DEC) R"""(>`` representation of pi. They can be compared afterwards to see if precision conversion is handled correctly.

:return: ``(string(pi) , boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<YADE_REAL_DEC>>(pi))`` from C++
	)""");

	py::def("testPiMpfr500", yade::testPi<yade::bmp::number<yade::bmp::mpfr_float_backend<500>>>, R"""(
This function returns a string and ``boost::mpfr_float_backend<500>`` representation of pi. They can be compared afterwards to see if precision conversion is handled correctly.

:return: ``(string(pi) , boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<500>>(pi))`` from C++
	)""");

#endif

}

