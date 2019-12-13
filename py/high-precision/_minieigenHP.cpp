/*************************************************************************
*  2009-2012 © Václav Šmilauer                                           *
*  2019        Janek Kozicki                                             *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "Real/Real.hpp"

#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

using namespace ::yade::MathEigenTypes;
#include <minieigen/converters.hpp>
#include <minieigen/visitors.hpp>
#include <minieigen/expose.hpp>

//#define ARBITRARY_REAL_DEBUG
#include "ExposeStorageOrdering.hpp"
#include "Real/ToFromPythonConverter.hpp"

BOOST_PYTHON_MODULE(THE_CPP_NAME)
try {
	// arbitrary Real specific stuff: start
	ArbitraryReal_from_python<Real>();
	py::to_python_converter<Real, ArbitraryReal_to_python<Real>>();
	expose_storage_ordering();
	// arbitrary Real specific stuff: end

	py::scope().attr("__doc__") = "miniEigen is wrapper for a small part of the `Eigen <http://eigen.tuxfamily.org>`_ library. Refer to its documentation "
	                              "for details. All classes in this module support pickling.";

	py::docstring_options docopt;
	docopt.enable_all();
	docopt.disable_cpp_signatures();


	expose_converters(); // in expose-converters.cpp

	expose_vectors();
	expose_matrices(); // must come after vectors
	expose_complex();
	expose_quaternion();
	expose_boxes();

	// Requires -ldouble-conversion, but that function isn't used anywhere
	//	py::def("float2str",&doubleToShortest,(py::arg("f"),py::arg("pad")=0),"Return the shortest string representation of *f* which will is equal to *f* when converted back to float. This function is only useful in Python prior to 3.0; starting from that version, standard string conversion does just that.");

#ifdef EIGEN_DONT_ALIGN
	py::scope().attr("vectorize") = false;
#else
	py::scope().attr("vectorize") = true;
#endif
} catch (...) {
	std::cerr << "Importing this module caused an exception and this module is in an inconsistent state now.\n\n";
	PyErr_Print();
	PyErr_SetString(PyExc_SystemError, __FILE__);
	boost::python::handle_exception();
	throw;
}

