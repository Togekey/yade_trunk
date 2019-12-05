/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// How I was checking std::complex:
//    https://cpp.hotexamples.com/examples/-/-/PyFloat_AsDouble/cpp-pyfloat_asdouble-function-examples.html
//    https://docs.python.org/3/c-api/float.html
//    https://bugs.python.org/issue628842
// How to use PyErr_Occurred()
//    https://stackoverflow.com/questions/1796510/accessing-a-python-traceback-from-the-c-api
//    https://docs.python.org/3/c-api/exceptions.html#c.PyErr_Occurred
// minieigen sources
//    apt-get source minieigen
//    /tmp/minieigen/minieigen-0.50.3+dfsg1/src
// How to call mpmath, tzn boost::python::eval
//    https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/tutorial/tutorial/embedding.html
//    https://docs.python.org/2/extending/index.html
//    https://www.boost.org/doc/libs/1_47_0/libs/python/doc/v2/import.html
//    https://stackoverflow.com/questions/937884/how-do-i-import-modules-in-boostpython-embedded-python-code
// !→ https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
//    https://docs.python.org/2/extending/embedding.html
// How I was grepping python sources:
//    /home/deb/deb/python/python3.7-3.7.3  z apt-get source (ctags)
// !→ /home/deb/deb/python/pybind11-2.2.4/include/pybind11  też z apt-get source (ctags)
//    https://packages.debian.org/search?suite=sid&searchon=names&keywords=pybind
//    https://github.com/pybind/pybind11
// Interesting mppp library:
//    https://github.com/bluescarni/mppp
//    https://bluescarni.github.io/mppp/
//    git clone https://github.com/bluescarni/mppp.git
//    /home/deb/deb/python/mppp/mppp   (ctags)
// !→ /home/deb/deb/python/mppp/mppp/include/mp++/extra/pybind11.hpp
// More links:
//    http://mpmath.org/    →     https://github.com/fredrik-johansson/mpmath/blob/master/mpmath/ctx_mp_python.py
//    http://flintlib.org/ który używa Arb i z niego korzysta mppp
//    https://packages.debian.org/search?searchon=sourcenames&keywords=flint
//    http://fredrikj.net/math/x2019.pdf  Arb
// !→ https://pyformat.info/
// Boost python:
//    https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/tutorial/index.html
//    https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/article.html
//    https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/index.html
//    https://www.boost.org/doc/libs/1_42_0/libs/python/doc/v2/faq.html#topythonconversionfailed
//    https://docs.python.org/2.5/ext/dnt-basics.html
//    https://wiki.python.org/moin/boost.python/HowTo

#include "Real/Real.hpp"

#include <boost/python.hpp>
#include <iostream>
#include <limits>
#include <sstream>

#define ARBITRARY_REAL_DEBUG
#include "Real/ToFromPythonConverter.hpp"

namespace py = ::boost::python;

struct Var {
	Real value { -71.23 };
	Real get() const { return value; };
	void set(Real val) { value = val; };
};

#include <boost/python/def.hpp>
int f(double x, double y, double z = 0.0, double w = 1.0)
{
	std::cerr << "x=" << x << "\n";
	std::cerr << "y=" << y << "\n";
	std::cerr << "z=" << z << "\n";
	std::cerr << "w=" << w << "\n";
	return 92;
}

double someFunction()
{
	double a   = 10.11;
	double b   = 30.11;
	double y   = a * b + b * a + b - a;
	double ret = std::pow(y, a);
	std::cerr << ret << "\n";
	return ret;
}

// ALL DOCS         : https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/index.html
//              note: clicking 'next →' skips tutorial and references.
//                    https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/tutorial
//                    https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/index.html
//              good: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/tutorial/tutorial/techniques.html
//                    https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/high_level_components.html
//         operators: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/high_level_components/boost_python_operators_hpp.html
//              enum: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/tutorial/tutorial/object.html#tutorial.object.enums
//                    https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/high_level_components/boost_python_enum_hpp.html
//             scope: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/high_level_components/boost_python_scope_hpp.html
// python object ref: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/object_wrappers/boost_python_object_hpp.html
//               arg: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/function_invocation_and_creation.html
// arbitrary arg num: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/function_invocation_and_creation/boost_python_raw_function_hpp.html
//  call python func: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/topics.html
//                    https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/function_invocation_and_creation/boost_python_call_hpp.html
//                    https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/function_invocation_and_creation/boost_python_call_method_hpp.html
// virtual functions: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/topics.html#topics.calling_python_functions_and_met.introduction
//               ptr: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/function_invocation_and_creation/boost_python_ptr_hpp.html
//                    https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/to_from_python_type_conversion/boost_python_register_ptr_to_pyt.html
//                    https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/utility_and_infrastructure/boost_python_pointee_hpp.html
//  handle, borrowed: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/utility_and_infrastructure/boost_python_handle_hpp.html
//            pickle: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/topics/pickle_support.html
//                    https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/topics/pickle_support.html#topics.pickle_support.light_weight_alternative_pickle_
//   indexing vector: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/topics/indexing_support.html
//       boost numpy: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/numpy/index.html
// docstring_options: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/function_invocation_and_creation/function_documentation.html
//                    https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/high_level_components/boost_python_docstring_options_h.html

BOOST_PYTHON_MODULE(THE_CPP_NAME)
try {
	// what about BOOST_PYTHON_PY_SIGNATURES_PROPER_INIT_SELF_TYPE
	// If defined the python type of __init__ method "self" parameters is properly generated, otherwise object is used.
	// https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/configuration.html

	// If it can't find a propeer conversion:
	// https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/faq/why_is_my_automatic_to_python_co.html

	ArbitraryReal_from_python<Real>();
	py::to_python_converter<Real, ArbitraryReal_to_python<Real>>();
#ifdef YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this
	mpfr::mpreal::set_default_prec(YADE_REAL_BIT + 1);
#endif

	Real startPrecision = std::numeric_limits<Real>::epsilon(); //Eigen::NumTraits<Scalar>::dummy_precision();
	std::cerr << __FILE__ << ":" << __LINE__ << " startPrecision: " << startPrecision << "\n";
	boost::python::object startPrecisionPyObj { startPrecision };

	py::class_<Var>("Var").add_property("val", &Var::get, &Var::set);

	py::def("f", f, (py::arg("x"), "y", py::arg("z") = 0.0, py::arg("w") = someFunction()));
	/*} catch(const py::error_already_set& e) {
	PyErr_Print();
//	boost::python::handle_exception();
throw e;
//	PyErr_Clear();
}
catch (...) {
	PyObject *type, *value, *traceback;
	PyErr_Fetch(&type, &value, &traceback);
	throw std::runtime_error(boost::python::call_method<std::string>(value, "__str__"));
	//PyErr_Print();
//	PyErr_Restore(type, value, traceback);
	//PyErr_SetImportError( boost::python::incref(boost::python::object(__FILE__).ptr()) , nullptr, nullptr);
	//boost::python::throw_error_already_set();
//	throw boost::python::error_already_set();
//	boost::python::handle_exception();
//	std::cerr << "Unknown exception\n";
//	throw;
}
*/
} catch (...) {
	// see: https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/high_level_components/boost_python_errors_hpp.html#high_level_components.boost_python_errors_hpp.example
	//      https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/tutorial/tutorial/embedding.html
	//      https://stackoverflow.com/questions/1418015/how-to-get-python-exception-text
	// If we wanted custom yade exceptions thrown to python:
	//     https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/tutorial/tutorial/exception.html
	//     https://www.boost.org/doc/libs/1_71_0/libs/python/doc/html/reference/high_level_components/boost_python_exception_translato.html
	std::cerr << ("Importing this module caused an unrecognized exception caught on C++ side and this module is in an inconsistent state now.\n\n");
	PyErr_Print();
	PyErr_SetString(PyExc_SystemError, __FILE__); // raising anything other than SystemError is not possible
	boost::python::handle_exception();
	throw;
}

