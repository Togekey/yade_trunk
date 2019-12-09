/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "Real/Real.hpp"

#include <boost/python.hpp>
#include <iostream>
#include <limits>
#include <sstream>
#include <Eigen/Core>

//#define ARBITRARY_REAL_DEBUG
#include "Real/ToFromPythonConverter.hpp"
#include "ExposeStorageOrdering.hpp"

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
	//	std::cerr << ret << " !!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	return ret;
}

#if not (defined(EIGEN_NUM_TRAITS_HPP) or defined(EIGEN_MPREALSUPPORT_MODULE_H) or defined(YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this))
namespace boost {
namespace multiprecision {
}
}

namespace Substitute {
static constexpr long get_default_prec = std::numeric_limits<Real>::digits;
static inline Real    highest(long = get_default_prec) { return std::numeric_limits<Real>::max(); }
static inline Real    lowest(long = get_default_prec) { return std::numeric_limits<Real>::lowest(); }
static inline Real    Pi(long = get_default_prec) { return boost::math::constants::pi<Real>(); }
static inline Real    Euler(long = get_default_prec) { return boost::math::constants::euler<Real>(); }
static inline Real    Log2(long = get_default_prec)
{
	using namespace boost::multiprecision;
	using namespace std;
	return log(Real(2));
}
static inline Real Catalan(long = get_default_prec) { return boost::math::constants::catalan<Real>(); }

static inline Real epsilon(long = get_default_prec) { return std::numeric_limits<Real>::epsilon(); }
static inline Real epsilon(const Real&) { return std::numeric_limits<Real>::epsilon(); }
inline bool        isEqualFuzzy(const Real& a, const Real& b, const Real& eps)
{
	using namespace boost::multiprecision;
	using namespace std;
	return abs(a - b) <= eps;
}
}
#endif

BOOST_PYTHON_MODULE(THE_CPP_NAME)
try {
#ifdef YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this
	mpfr::mpreal::set_default_prec(YADE_REAL_BIT + 1);
	long defprec  = mpfr::mpreal::get_default_prec();
	long max_exp2 = mpfr::mpreal::get_emax();
#else
	long defprec  = std::numeric_limits<Real>::digits;
	long max_exp2 = std::numeric_limits<Real>::max_exponent;
#endif
	ArbitraryReal_from_python<Real>();
	py::to_python_converter<Real, ArbitraryReal_to_python<Real>>();

	py::class_<Var>("Var").add_property("val", &Var::get, &Var::set);

	py::def("f", f, (py::arg("x"), "y", py::arg("z") = 0.0, py::arg("w") = someFunction()));

	py::scope().attr("defprec")  = defprec;
	py::scope().attr("max_exp2") = max_exp2;

	// by zadziałało, gdyby ten enum miał nazwę.
	//py::enum_<Eigen::NumTraits<Real>::NazwaEnuma>("traits").value("IsInteger",Eigen::NumTraits<Real>::IsInteger).export_values();

	py::scope().attr("IsInteger")             = int(Eigen::NumTraits<Real>::IsInteger);
	py::scope().attr("IsSigned")              = int(Eigen::NumTraits<Real>::IsSigned);
	py::scope().attr("IsComplex")             = int(Eigen::NumTraits<Real>::IsComplex);
	py::scope().attr("RequireInitialization") = int(Eigen::NumTraits<Real>::RequireInitialization);
	py::scope().attr("ReadCost")              = int(Eigen::NumTraits<Real>::ReadCost);
	py::scope().attr("AddCost")               = int(Eigen::NumTraits<Real>::AddCost);
	py::scope().attr("MulCost")               = int(Eigen::NumTraits<Real>::MulCost);

#if defined(EIGEN_NUM_TRAITS_HPP) or defined(EIGEN_MPREALSUPPORT_MODULE_H) or defined(YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this)
	py::def("highest", Eigen::NumTraits<Real>::highest, (py::arg("Precision") = defprec));
	py::def("lowest", Eigen::NumTraits<Real>::lowest, (py::arg("Precision") = defprec));

	py::def("Pi", Eigen::NumTraits<Real>::Pi, (py::arg("Precision") = defprec));
	py::def("Euler", Eigen::NumTraits<Real>::Euler, (py::arg("Precision") = defprec));
	py::def("Log2", Eigen::NumTraits<Real>::Log2, (py::arg("Precision") = defprec));
	py::def("Catalan", Eigen::NumTraits<Real>::Catalan, (py::arg("Precision") = defprec));

	py::def("epsilon", static_cast<Real (*)(long)>(&Eigen::NumTraits<Real>::epsilon), (py::arg("Precision") = defprec));
	py::def("epsilon", static_cast<Real (*)(const Real&)>(&Eigen::NumTraits<Real>::epsilon), (py::arg("x")));
#ifdef YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this
	py::def("isEqualFuzzy", static_cast<bool (*)(const Real&, const Real&, const Real&)>(&mpfr::isEqualFuzzy));
#else
	py::def("isEqualFuzzy", Eigen::internal::isEqualFuzzy);
#endif
#else
	py::def("highest", Substitute::highest, (py::arg("Precision") = defprec));
	py::def("lowest", Substitute::lowest, (py::arg("Precision") = defprec));

	py::def("Pi", Substitute::Pi, (py::arg("Precision") = defprec));
	py::def("Euler", Substitute::Euler, (py::arg("Precision") = defprec));
	py::def("Log2", Substitute::Log2, (py::arg("Precision") = defprec));
	py::def("Catalan", Substitute::Catalan, (py::arg("Precision") = defprec));

	py::def("epsilon", static_cast<Real (*)(long)>(&Substitute::epsilon), (py::arg("Precision") = defprec));
	py::def("epsilon", static_cast<Real (*)(const Real&)>(&Substitute::epsilon), (py::arg("x")));
	py::def("isEqualFuzzy", Substitute::isEqualFuzzy);
#endif

	py::def("dummy_precision", Eigen::NumTraits<Real>::dummy_precision);

	py::def("random", static_cast<Real (*)()>(&Eigen::internal::random<Real>));
	py::def("random", static_cast<Real (*)(const Real&, const Real&)>(&Eigen::internal::random<Real>), (py::arg("a"), "b"));
	py::def("isMuchSmallerThan",
	        static_cast<bool (*)(const Real&, const Real&, const Real&)>(&Eigen::internal::isMuchSmallerThan),
	        (py::arg("a"), "b", "eps"));
	py::def("isApprox", static_cast<bool (*)(const Real&, const Real&, const Real&)>(&Eigen::internal::isApprox), (py::arg("a"), "b", "eps"));
	py::def("isApproxOrLessThan",
	        static_cast<bool (*)(const Real&, const Real&, const Real&)>(&Eigen::internal::isApproxOrLessThan),
	        (py::arg("a"), "b", "eps"));

	py::def("toLongDouble", static_cast<long double (*)(const Real&)>(&Eigen::internal::cast<Real, long double>), (py::arg("x")));
	py::def("toDouble", static_cast<double (*)(const Real&)>(&Eigen::internal::cast<Real, double>), (py::arg("x")));
	py::def("toLong", static_cast<long (*)(const Real&)>(&Eigen::internal::cast<Real, long>), (py::arg("x")));
	py::def("toInt", static_cast<int (*)(const Real&)>(&Eigen::internal::cast<Real, int>), (py::arg("x")));

	expose_storage_ordering();
} catch (...) {
	std::cerr << ("Importing this module caused an unrecognized exception caught on C++ side and this module is in an inconsistent state now.\n\n");
	PyErr_Print();
	PyErr_SetString(PyExc_SystemError, __FILE__); // raising anything other than SystemError is not possible
	boost::python::handle_exception();
	throw;
}

