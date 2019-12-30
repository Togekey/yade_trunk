/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// This python module exposes all C++ math functions for Real and Complex type to python.
// In fact it "just duplicates" 'import math', 'import cmath' or 'import mpmath'.
// This module has three purposes:
// 1. to reliably test all C++ math functions of arbitrary Real and Complex types against mpmath.
// 2. to test Eigen NumTraits
// 3. To allow writing python math code in  a way that mirrors C++ math code in yade. As a bonus it will be faster than mpmath
//    because mpmath is a purely python library (which was one of the main difficulties when writing lib/high-precision/ToFromPythonConverter.hpp)

// TODO: add docstrings to all these math functions.

#include <lib/high-precision/Real.hpp>
#ifdef YADE_CGAL
#include <lib/base/AliasCGAL.hpp>
#endif

#include <Eigen/Core>
#include <Eigen/src/Core/MathFunctions.h>
#include <boost/python.hpp>
#include <iostream>
#include <limits>
#include <sstream>

#include <py/high-precision/_ExposeStorageOrdering.hpp>

// testing Real type
#include <boost/concept/assert.hpp>
#include <boost/math/concepts/real_type_concept.hpp>

namespace py = ::boost::python;
using ::yade::Complex;
using ::yade::Real;

// Converts a std::pair instance to a Python tuple.
template <typename T1, typename T2> struct std_pair_to_tuple {
	static PyObject*           convert(std::pair<T1, T2> const& p) { return boost::python::incref(boost::python::make_tuple(p.first, p.second).ptr()); }
	static PyTypeObject const* get_pytype() { return &PyTuple_Type; }
};

// Helper for convenience.
template <typename T1, typename T2> struct std_pair_to_python_converter {
	std_pair_to_python_converter()
	{
		boost::python::to_python_converter<
		        std::pair<T1, T2>,
		        std_pair_to_tuple<T1, T2>,
		        true //std_pair_to_tuple has get_pytype
		        >();
	}
};

std::pair<Real, int> test_frexp(const Real& x)
{
	int  i   = 0;
	Real ret = ::yade::math::frexp(x, &i);
	return std::pair<Real, int> { ret, i };
}

std::pair<Real, Real> test_modf(const Real& x)
{
	Real r = 0;
#if defined(YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this)
	Real ret = ::yade::math::modf(x, r);
#else
	Real ret      = ::yade::math::modf(x, &r);
#endif
	return std::pair<Real, Real> { ret, r };
}

std::pair<Real, long> test_remquo(const Real& x, const Real& y)
{
#if defined(YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this)
	long i = 0;
#else
	int  i        = 0;
#endif
	Real ret = ::yade::math::remquo(x, y, &i);
	return std::pair<Real, long> { ret, i };
}


#ifdef YADE_CGAL

bool test_CGAL_Is_valid(const Real& x) {
	return CGAL::Is_valid<Real>()(x);
}

Real test_CGAL_Sqrt(const Real& x) {
	return CGAL::Algebraic_structure_traits<Real>::Sqrt()(x);
}

Real test_CGAL_Kth_root(int k, const Real& x) {
	return CGAL::Algebraic_structure_traits<Real>::Kth_root()(k,x);
}

std::pair<double, double> test_CGAL_To_interval(const Real& x) {
	return CGAL::Real_embeddable_traits<Real>::To_interval()(x);
}

bool test_CGAL_Is_finite(const Real& x) {
	return CGAL::Real_embeddable_traits<Real>::Is_finite()(x);
}

#endif



struct Var {
	Real    value { -71.23 };
	Complex valueComplex { -71.23, 33.23 };

	Real get() const { return value; };
	void set(Real val) { value = val; };

	Complex getComplex() const { return valueComplex; };
	void    setComplex(Complex val) { valueComplex = val; };
};

void compareVec(const std::vector<Real>& vec, const ::yade::math::UnderlyingReal* array)
{
	for (size_t i = 0; i < vec.size(); i++) {
		if (vec[i] != array[i]) {
			std::cerr << __PRETTY_FUNCTION__ << " failed test\n";
			exit(1);
		}
	}
}

#include <boost/range/combine.hpp>
// this funcction simulates some external library which works on C-arrays.
void multVec(::yade::math::UnderlyingReal* array, const ::yade::math::UnderlyingReal& fac, size_t s)
{
	for (size_t i = 0; i < s; i++)
		array[i] *= fac;
}

namespace yade {
void testArray()
{
	std::vector<Real> vec {};
	int               i = 1000;
	while (i-- > 0)
		vec.push_back(math::random01());
	compareVec(vec, math::constVectorData(vec));
	auto copy = vec;
	Real fac  = 3.33;
	multVec(math::vectorData(vec), static_cast<::yade::math::UnderlyingReal>(fac), vec.size());
	for (auto a : boost::combine(copy, vec)) {
		if (a.get<0>() * fac != a.get<1>()) {
			std::cerr << __PRETTY_FUNCTION__ << " failed test\n";
			exit(1);
		}
	}
}
}

namespace {
static inline Real smallest_positive() { return std::numeric_limits<Real>::min(); }
}

#if not(defined(YADE_EIGEN_NUM_TRAITS_HPP) or defined(EIGEN_MPREALSUPPORT_MODULE_H) or defined(YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this))
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

BOOST_PYTHON_MODULE(_math)
try {
#ifdef YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this
	mpfr::mpreal::set_default_prec(YADE_REAL_BIT + 1);
	long defprec  = mpfr::mpreal::get_default_prec();
	long max_exp2 = mpfr::mpreal::get_emax();
#else
	long defprec  = std::numeric_limits<Real>::digits;
	long max_exp2 = std::numeric_limits<Real>::max_exponent;
#endif
	// This is registered in py/high-precision/_minieigenHP.cpp
	//ArbitraryComplex_from_python<Complex>();
	//py::to_python_converter<Complex, ArbitraryComplex_to_python<Complex>>();

	//ArbitraryReal_from_python<Real>();
	//py::to_python_converter<Real, ArbitraryReal_to_python<Real>>();

	py::class_<Var>("Var").add_property("val", &Var::get, &Var::set).add_property("cpl", &Var::getComplex, &Var::setComplex, "The ``Var`` class is used to test to/from python converters for arbitrary precision ``Real``. It has one ``Real`` (``val``) and one ``Complex`` (``cpl``) variable to test reading from and writing to.");

	py::scope().attr("defprec")  = defprec;
	py::scope().attr("max_exp2") = max_exp2;

	// it would work if this enum had a name.
	//py::enum_<Eigen::NumTraits<Real>::NazwaEnuma>("traits").value("IsInteger",Eigen::NumTraits<Real>::IsInteger).export_values();

	py::scope().attr("IsInteger")             = int(Eigen::NumTraits<Real>::IsInteger);
	py::scope().attr("IsSigned")              = int(Eigen::NumTraits<Real>::IsSigned);
	py::scope().attr("IsComplex")             = int(Eigen::NumTraits<Real>::IsComplex);
	py::scope().attr("RequireInitialization") = int(Eigen::NumTraits<Real>::RequireInitialization);
	py::scope().attr("ReadCost")              = int(Eigen::NumTraits<Real>::ReadCost);
	py::scope().attr("AddCost")               = int(Eigen::NumTraits<Real>::AddCost);
	py::scope().attr("MulCost")               = int(Eigen::NumTraits<Real>::MulCost);

#if defined(YADE_EIGEN_NUM_TRAITS_HPP) or defined(EIGEN_MPREALSUPPORT_MODULE_H) or defined(YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this)
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
	py::def("smallest_positive", static_cast<Real (*)()>(&smallest_positive));
#else
	py::def("isEqualFuzzy", Eigen::internal::isEqualFuzzy);
	py::def("smallest_positive", static_cast<Real (*)()>(&Eigen::NumTraits<Real>::smallest_positive));
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
	py::def("smallest_positive", static_cast<Real (*)()>(&smallest_positive));
	py::def("isEqualFuzzy", Substitute::isEqualFuzzy);
#endif

	py::def("dummy_precision", Eigen::NumTraits<Real>::dummy_precision);

	py::def("random", static_cast<Real (*)()>(&Eigen::internal::random<Real>));
	py::def("random", static_cast<Real (*)(const Real&, const Real&)>(&Eigen::internal::random<Real>), (py::arg("a"), "b"));

#if ((EIGEN_MAJOR_VERSION > 2) and (EIGEN_WORLD_VERSION >= 3)) or defined(YADE_EIGEN_NUM_TRAITS_HPP) or defined(EIGEN_MPREALSUPPORT_MODULE_H)                  \
        or defined(YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this)
	py::def("isMuchSmallerThan",
	        static_cast<bool (*)(const Real&, const Real&, const Real&)>(&Eigen::internal::isMuchSmallerThan),
	        (py::arg("a"), "b", "eps"));
	py::def("isApprox", static_cast<bool (*)(const Real&, const Real&, const Real&)>(&Eigen::internal::isApprox), (py::arg("a"), "b", "eps"));
	py::def("isApproxOrLessThan",
	        static_cast<bool (*)(const Real&, const Real&, const Real&)>(&Eigen::internal::isApproxOrLessThan),
	        (py::arg("a"), "b", "eps"));
#else
	// older eigen 3.2 didn't use `const Real&` but was copying third argument by value `Real`
	py::def("isMuchSmallerThan",
	        static_cast<bool (*)(const Real&, const Real&, Real)>(&Eigen::internal::isMuchSmallerThan<Real, Real>),
	        (py::arg("a"), "b", "eps"));
	py::def("isApprox", static_cast<bool (*)(const Real&, const Real&, Real)>(&Eigen::internal::isApprox<Real>), (py::arg("a"), "b", "eps"));
	py::def("isApproxOrLessThan",
	        static_cast<bool (*)(const Real&, const Real&, Real)>(&Eigen::internal::isApproxOrLessThan<Real>),
	        (py::arg("a"), "b", "eps"));
#endif

	py::def("toLongDouble", static_cast<long double (*)(const Real&)>(&Eigen::internal::cast<Real, long double>), (py::arg("x")));
	py::def("toDouble", static_cast<double (*)(const Real&)>(&Eigen::internal::cast<Real, double>), (py::arg("x")));
	py::def("toLong", static_cast<long (*)(const Real&)>(&Eigen::internal::cast<Real, long>), (py::arg("x")));
	py::def("toInt", static_cast<int (*)(const Real&)>(&Eigen::internal::cast<Real, int>), (py::arg("x")));

	expose_storage_ordering();

#ifdef YADE_REAL_MPFR_NO_BOOST_experiments_only_never_use_this
#warning "::mpfr::mpreal (non-boost implementation) is not passing Boost RealTypeConcept test"
#else
	BOOST_CONCEPT_ASSERT((boost::math::concepts::RealTypeConcept<Real>));
#endif

	// check overload (and namespace) resolution for all math functions. As a side effect they are exported to python, and can be unit-tested.
#define YADE_PYEXPORT_MATH_1(func) py::def(#func, static_cast<Real (*)(const Real&)>(&::yade::math::func), (py::arg("x")));
#define YADE_PYEXPORT_MATH_1_MREF(func) py::def(#func, static_cast<Real (*)(Real)>(&::yade::math::func), (py::arg("x")));
#define YADE_PYEXPORT_MATH_1_COMPLEX(func) py::def(#func, static_cast<Complex (*)(const Complex&)>(&::yade::math::func), (py::arg("x")));
#define YADE_PYEXPORT_MATH_1_COMPLEX_TO_REAL(func) py::def(#func, static_cast<Real (*)(const Complex&)>(&::yade::math::func), (py::arg("x")));
#define YADE_PYEXPORT_MATH_1_INT(func) py::def(#func, static_cast<int (*)(const Real&)>(&::yade::math::func), (py::arg("x")));
	// FIXED: maybe registering the complex versions first will solve the problem that they are used when Real ones should be used - Yes.
	// Complex versions are registered first
	YADE_PYEXPORT_MATH_1_COMPLEX(sin)
	YADE_PYEXPORT_MATH_1_COMPLEX(sinh)
	YADE_PYEXPORT_MATH_1_COMPLEX(cos)
	YADE_PYEXPORT_MATH_1_COMPLEX(cosh)
	YADE_PYEXPORT_MATH_1_COMPLEX(tan)
	YADE_PYEXPORT_MATH_1_COMPLEX(tanh)

	YADE_PYEXPORT_MATH_1_COMPLEX_TO_REAL(abs)
	YADE_PYEXPORT_MATH_1_COMPLEX(conj)
	YADE_PYEXPORT_MATH_1_COMPLEX_TO_REAL(real)
	YADE_PYEXPORT_MATH_1_COMPLEX_TO_REAL(imag)

	YADE_PYEXPORT_MATH_1_COMPLEX(exp)
	YADE_PYEXPORT_MATH_1_COMPLEX(log)

	// Real versions are registered afterwards
	YADE_PYEXPORT_MATH_1(sin)
	YADE_PYEXPORT_MATH_1(sinh)
	YADE_PYEXPORT_MATH_1(cos)
	YADE_PYEXPORT_MATH_1(cosh)
	YADE_PYEXPORT_MATH_1(tan)
	YADE_PYEXPORT_MATH_1(tanh)

#if (YADE_REAL_BIT > 64)
	YADE_PYEXPORT_MATH_1(abs)
	YADE_PYEXPORT_MATH_1(fabs)
#else
	YADE_PYEXPORT_MATH_1_MREF(fabs)
	YADE_PYEXPORT_MATH_1_MREF(abs)
#endif
	YADE_PYEXPORT_MATH_1(acos)
	YADE_PYEXPORT_MATH_1(acosh)
	YADE_PYEXPORT_MATH_1(asin)
	YADE_PYEXPORT_MATH_1(asinh)
	YADE_PYEXPORT_MATH_1(atan)
	YADE_PYEXPORT_MATH_1(atanh)
	YADE_PYEXPORT_MATH_1(cbrt)
	YADE_PYEXPORT_MATH_1(ceil)
	YADE_PYEXPORT_MATH_1(erf)
	YADE_PYEXPORT_MATH_1(erfc)
	YADE_PYEXPORT_MATH_1(exp)
	YADE_PYEXPORT_MATH_1(exp2)
	YADE_PYEXPORT_MATH_1(expm1)
	YADE_PYEXPORT_MATH_1(floor)
	YADE_PYEXPORT_MATH_1(ilogb)
	YADE_PYEXPORT_MATH_1(lgamma)
	YADE_PYEXPORT_MATH_1(log)
	YADE_PYEXPORT_MATH_1(log10)
	YADE_PYEXPORT_MATH_1(log1p)
	YADE_PYEXPORT_MATH_1(log2)
	YADE_PYEXPORT_MATH_1(logb)
	//YADE_PYEXPORT_MATH_1(riemann_zeta) // since C++17
	YADE_PYEXPORT_MATH_1(rint)
	YADE_PYEXPORT_MATH_1(round)
	YADE_PYEXPORT_MATH_1(sqrt)
	YADE_PYEXPORT_MATH_1(tgamma)
	YADE_PYEXPORT_MATH_1(trunc)

	YADE_PYEXPORT_MATH_1_INT(sgn)
	YADE_PYEXPORT_MATH_1_INT(sign)
#undef YADE_PYEXPORT_MATH_1
#undef YADE_PYEXPORT_MATH_1_MREF
#undef YADE_PYEXPORT_MATH_1_COMPLEX
#undef YADE_PYEXPORT_MATH_1_COMPLEX_TO_REAL
#undef YADE_PYEXPORT_MATH_1_INT

#define YADE_PYEXPORT_MATH_2(func) py::def(#func, static_cast<Real (*)(const Real&, const Real&)>(&::yade::math::func), (py::arg("x"), "y"));
#define YADE_PYEXPORT_MATH_2_CREF(func)                                                                                                                        \
	py::def(#func,                                                                                                                                         \
	        static_cast<const Real& (*)(const Real&, const Real&)>(&::yade::math::func),                                                                   \
	        (py::arg("x"), "y"),                                                                                                                           \
	        py::return_value_policy<py::copy_const_reference>());
	YADE_PYEXPORT_MATH_2(atan2)
	//YADE_PYEXPORT_MATH_2(beta) // since C++17
	//YADE_PYEXPORT_MATH_2(cyl_bessel_i) // since C++17
	//YADE_PYEXPORT_MATH_2(cyl_bessel_j) // since C++17
	//YADE_PYEXPORT_MATH_2(cyl_bessel_k) // since C++17
	YADE_PYEXPORT_MATH_2(fmod)
	YADE_PYEXPORT_MATH_2(hypot)
	YADE_PYEXPORT_MATH_2_CREF(max)
	YADE_PYEXPORT_MATH_2_CREF(min)
	YADE_PYEXPORT_MATH_2(pow)
	YADE_PYEXPORT_MATH_2(remainder)
#undef YADE_PYEXPORT_MATH_2
#undef YADE_PYEXPORT_MATH_2_CREF

#define YADE_PYEXPORT_MATH_2_TYPE1(func, FirstType) py::def(#func, static_cast<Real (*)(FirstType, const Real&)>(&::yade::math::func), (py::arg("x"), "y"));
	//YADE_PYEXPORT_MATH_2_TYPE1(sph_bessel, unsigned) // since C++17
#undef YADE_PYEXPORT_MATH_2_TYPE1

#define YADE_PYEXPORT_MATH_2_TYPE2(func, SecondType) py::def(#func, static_cast<Real (*)(const Real&, SecondType)>(&::yade::math::func), (py::arg("x"), "y"));
	YADE_PYEXPORT_MATH_2_TYPE2(ldexp, int)
#undef YADE_PYEXPORT_MATH_2_TYPE2

	std_pair_to_python_converter<Real, Real>();
	std_pair_to_python_converter<double, double>();
	std_pair_to_python_converter<Real, long>();
	std_pair_to_python_converter<Real, int>();
	py::def("frexp", test_frexp, (py::arg("x")));
	py::def("modf", test_modf, (py::arg("x")));

#define YADE_PYEXPORT_MATH_3(func) py::def(#func, static_cast<Real (*)(const Real&, const Real&, const Real&)>(&::yade::math::func), (py::arg("x"), "y", "z"));
	YADE_PYEXPORT_MATH_3(fma)
#undef YADE_PYEXPORT_MATH_3

	py::def("remquo", test_remquo, (py::arg("x"), "y"));
	py::def("testArray", ::yade::testArray);

#ifdef YADE_CGAL
	py::scope().attr("testCgalNumTraits") = true;
	py::def("CGAL_Is_valid", test_CGAL_Is_valid, (py::arg("x")));
	py::def("CGAL_Sqrt", test_CGAL_Sqrt, (py::arg("x")));
	py::def("CGAL_Kth_root", test_CGAL_Kth_root, (py::arg("x")));
	py::def("CGAL_To_interval", test_CGAL_To_interval, (py::arg("x")));
	py::def("CGAL_Is_finite", test_CGAL_Is_finite, (py::arg("x")));
#else
	py::scope().attr("testCgalNumTraits") = false;
#endif

} catch (...) {
	std::cerr << ("Importing this module caused an unrecognized exception caught on C++ side and this module is in an inconsistent state now.\n\n");
	PyErr_Print();
	PyErr_SetString(PyExc_SystemError, __FILE__); // raising anything other than SystemError is not possible
	boost::python::handle_exception();
	throw;
}

