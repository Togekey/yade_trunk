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

#include <lib/base/Logging.hpp>
#include <lib/high-precision/Real.hpp>
#include <lib/pyutil/doc_opts.hpp>
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

CREATE_CPP_LOCAL_LOGGER("_math.cpp")

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
	Real r   = 0;
	Real ret = ::yade::math::modf(x, &r);
	return std::pair<Real, Real> { ret, r };
}

std::pair<Real, long> test_remquo(const Real& x, const Real& y)
{
	int  i   = 0;
	Real ret = ::yade::math::remquo(x, y, &i);
	return std::pair<Real, long> { ret, i };
}


#ifdef YADE_CGAL

bool                      test_CGAL_Is_valid(const Real& x) { return CGAL::Is_valid<Real>()(x); }
Real                      test_CGAL_Square(const Real& x) { return CGAL::Algebraic_structure_traits<Real>::Square()(x); }
Real                      test_CGAL_Sqrt(const Real& x) { return CGAL::Algebraic_structure_traits<Real>::Sqrt()(x); }
Real                      test_CGAL_Kth_root(int k, const Real& x) { return CGAL::Algebraic_structure_traits<Real>::Kth_root()(k, x); }
std::pair<double, double> test_CGAL_To_interval(const Real& x) { return CGAL::Real_embeddable_traits<Real>::To_interval()(x); }
int                       test_CGAL_Sgn(const Real& x) { return int(CGAL::Real_embeddable_traits<Real>::Sgn()(x)); }
bool                      test_CGAL_Is_finite(const Real& x) { return CGAL::Real_embeddable_traits<Real>::Is_finite()(x); }


namespace yade {
Real testSimpleCgalNumTraitsCalculation()
{
	CGALpoint  x(Real(1), Real(1), Real(1));
	CGALpoint  p1(Real(0), Real(0), Real(0));
	CGALvector v1(Real(1), Real(1), Real(1));
	Plane      P(p1, v1);
	Real h = P.a() * x.x() + P.b() * x.y() + P.c() * x.z() + P.d();
	return ((h > 0.) - (h < 0.)) * pow(h, 2) / (CGALvector(P.a(), P.b(), P.c())).squared_length();
}
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

#if not(defined(YADE_EIGEN_NUM_TRAITS_HPP) or defined(EIGEN_MPREALSUPPORT_MODULE_H))
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
	YADE_SET_DOCSTRING_OPTS;

	long defprec  = std::numeric_limits<Real>::digits;
	long max_exp2 = std::numeric_limits<Real>::max_exponent;
	// This is registered in py/high-precision/_minieigenHP.cpp
	//ArbitraryComplex_from_python<Complex>();
	//py::to_python_converter<Complex, ArbitraryComplex_to_python<Complex>>();

	//ArbitraryReal_from_python<Real>();
	//py::to_python_converter<Real, ArbitraryReal_to_python<Real>>();

	py::class_<Var>("Var", "The ``Var`` class is used to test to/from python converters for arbitrary precision ``Real``")
	        .add_property("val", &Var::get, &Var::set, "one ``Real`` variable for testing.")
	        .add_property("cpl", &Var::getComplex, &Var::setComplex, "one ``Complex`` variable to test reading from and writing to it.");

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

#if defined(YADE_EIGEN_NUM_TRAITS_HPP) or defined(EIGEN_MPREALSUPPORT_MODULE_H)
	py::def("highest",
	        Eigen::NumTraits<Real>::highest,
	        (py::arg("Precision") = defprec),
	        R"""(:return: ``Real`` returns the largest finite value of the ``Real`` type. Wraps `std::numeric_limits<Real>::max() <https://en.cppreference.com/w/cpp/types/numeric_limits/max>`__ function.)""");
	py::def("lowest",
	        Eigen::NumTraits<Real>::lowest,
	        (py::arg("Precision") = defprec),
	        R"""(:return: ``Real`` returns the lowest (negative) finite value of the ``Real`` type. Wraps `std::numeric_limits<Real>::lowest() <https://en.cppreference.com/w/cpp/types/numeric_limits/lowest>`__ function.)""");

	py::def("Pi",
	        Eigen::NumTraits<Real>::Pi,
	        (py::arg("Precision") = defprec),
	        R"""(:return: ``Real`` The `π constant <https://en.wikipedia.org/wiki/Pi>`__, exposed to python for :ysrc:`testing <py/tests/testMath.py>` of :ysrc:`eigen numerical traits<lib/high-precision/EigenNumTraits.hpp>`.)""");
	py::def("Euler",
	        Eigen::NumTraits<Real>::Euler,
	        (py::arg("Precision") = defprec),
	        R"""(:return: ``Real`` The `Euler–Mascheroni constant <https://en.wikipedia.org/wiki/Euler%E2%80%93Mascheroni_constant>`__, exposed to python for :ysrc:`testing <py/tests/testMath.py>` of :ysrc:`eigen numerical traits<lib/high-precision/EigenNumTraits.hpp>`.)""");
	py::def("Log2",
	        Eigen::NumTraits<Real>::Log2,
	        (py::arg("Precision") = defprec),
	        R"""(:return: ``Real`` natural logarithm of 2, exposed to python for :ysrc:`testing <py/tests/testMath.py>` of :ysrc:`eigen numerical traits<lib/high-precision/EigenNumTraits.hpp>`.)""");
	py::def("Catalan",
	        Eigen::NumTraits<Real>::Catalan,
	        (py::arg("Precision") = defprec),
	        R"""(:return: ``Real`` The `catalan constant <https://en.wikipedia.org/wiki/Catalan%27s_constant>`__, exposed to python for :ysrc:`testing <py/tests/testMath.py>` of :ysrc:`eigen numerical traits<lib/high-precision/EigenNumTraits.hpp>`.)""");

	py::def("epsilon",
	        static_cast<Real (*)(long)>(&Eigen::NumTraits<Real>::epsilon),
	        (py::arg("Precision") = defprec),
	        R"""(:return: ``Real`` returns the difference between ``1.0`` and the next representable value of the ``Real`` type. Wraps `std::numeric_limits<Real>::epsilon() <https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon>`__ function.)""");
	py::def("epsilon",
	        static_cast<Real (*)(const Real&)>(&Eigen::NumTraits<Real>::epsilon),
	        (py::arg("x")),
	        R"""(:return: ``Real`` returns the difference between ``1.0`` and the next representable value of the ``Real`` type. Wraps `std::numeric_limits<Real>::epsilon() <https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon>`__ function.)""");
	py::def("isEqualFuzzy",
	        Eigen::internal::isEqualFuzzy,
	        R"""(:return: ``bool``, ``True`` if the absolute difference between two numbers is smaller than `std::numeric_limits<Real>::epsilon() <https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon>`__)""");
	py::def("smallest_positive",
	        static_cast<Real (*)()>(&Eigen::NumTraits<Real>::smallest_positive),
	        R"""(:return: ``Real`` the smallest number greater than zero. Wraps `std::numeric_limits<Real>::min() <https://en.cppreference.com/w/cpp/types/numeric_limits/min>`__)""");
#else
	py::def("highest",
	        Substitute::highest,
	        (py::arg("Precision") = defprec),
	        R"""(:return: ``Real`` returns the largest finite value of the ``Real`` type. Wraps `std::numeric_limits<Real>::max() <https://en.cppreference.com/w/cpp/types/numeric_limits/max>`__ function.)""");
	py::def("lowest",
	        Substitute::lowest,
	        (py::arg("Precision") = defprec),
	        R"""(:return: ``Real`` returns the lowest (negative) finite value of the ``Real`` type. Wraps `std::numeric_limits<Real>::lowest() <https://en.cppreference.com/w/cpp/types/numeric_limits/lowest>`__ function.)""");

	py::def("Pi",
	        Substitute::Pi,
	        (py::arg("Precision") = defprec),
	        R"""(:return: ``Real`` The `π constant <https://en.wikipedia.org/wiki/Pi>`__, exposed to python for :ysrc:`testing <py/tests/testMath.py>` of :ysrc:`eigen numerical traits<lib/high-precision/EigenNumTraits.hpp>`.)""");
	py::def("Euler",
	        Substitute::Euler,
	        (py::arg("Precision") = defprec),
	        R"""(:return: ``Real`` The `Euler–Mascheroni constant <https://en.wikipedia.org/wiki/Euler%E2%80%93Mascheroni_constant>`__, exposed to python for :ysrc:`testing <py/tests/testMath.py>` of :ysrc:`eigen numerical traits<lib/high-precision/EigenNumTraits.hpp>`.)""");
	py::def("Log2",
	        Substitute::Log2,
	        (py::arg("Precision") = defprec),
	        R"""(:return: ``Real`` natural logarithm of 2, exposed to python for :ysrc:`testing <py/tests/testMath.py>` of :ysrc:`eigen numerical traits<lib/high-precision/EigenNumTraits.hpp>`.)""");
	py::def("Catalan",
	        Substitute::Catalan,
	        (py::arg("Precision") = defprec),
	        R"""(:return: ``Real`` The `catalan constant <https://en.wikipedia.org/wiki/Catalan%27s_constant>`__, exposed to python for :ysrc:`testing <py/tests/testMath.py>` of :ysrc:`eigen numerical traits<lib/high-precision/EigenNumTraits.hpp>`.)""");

	py::def("epsilon",
	        static_cast<Real (*)(long)>(&Substitute::epsilon),
	        (py::arg("Precision") = defprec),
	        R"""(:return: ``Real`` returns the difference between ``1.0`` and the next representable value of the ``Real`` type. Wraps `std::numeric_limits<Real>::epsilon() <https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon>`__ function.)""");
	py::def("epsilon",
	        static_cast<Real (*)(const Real&)>(&Substitute::epsilon),
	        (py::arg("x")),
	        R"""(:return: ``Real`` returns the difference between ``1.0`` and the next representable value of the ``Real`` type. Wraps `std::numeric_limits<Real>::epsilon() <https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon>`__ function.)""");
	py::def("smallest_positive",
	        static_cast<Real (*)()>(&smallest_positive),
	        R"""(:return: ``Real`` the smallest number greater than zero. Wraps `std::numeric_limits<Real>::min() <https://en.cppreference.com/w/cpp/types/numeric_limits/min>`__)""");
	py::def("isEqualFuzzy",
	        Substitute::isEqualFuzzy,
	        R"""(:return: ``bool``, ``True`` if the absolute difference between two numbers is smaller than `std::numeric_limits<Real>::epsilon() <https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon>`__)""");
#endif

	py::def("dummy_precision",
	        Eigen::NumTraits<Real>::dummy_precision,
	        R"""(:return: similar to the function ``epsilon``, but assumes that last 10% of bits contain the numerical error only. This is sometimes used by Eigen when calling ``isEqualFuzzy`` to determine if values differ a lot or if they are vaguely close to each other.)""");

	py::def("random",
	        static_cast<Real (*)()>(&Eigen::internal::random<Real>),
	        R"""(:return: ``Real`` a symmetric random number in interval ``(-1,1)``. Used by Eigen.)""");
	py::def("random",
	        static_cast<Real (*)(const Real&, const Real&)>(&Eigen::internal::random<Real>),
	        (py::arg("a"), "b"),
	        R"""(:return: ``Real`` a random number in interval ``(a,b)``. Used by Eigen.)""");

#if ((EIGEN_MAJOR_VERSION > 2) and (EIGEN_WORLD_VERSION >= 3)) or defined(YADE_EIGEN_NUM_TRAITS_HPP) or defined(EIGEN_MPREALSUPPORT_MODULE_H)
	py::def("isMuchSmallerThan",
	        static_cast<bool (*)(const Real&, const Real&, const Real&)>(&Eigen::internal::isMuchSmallerThan),
	        (py::arg("a"), "b", "eps"),
	        R"""(:return: ``bool``, True if ``a`` is less than ``b`` with provided ``eps``, see also `here <https://stackoverflow.com/questions/15051367/how-to-compare-vectors-approximately-in-eigen>`__)""");
	py::def("isApprox",
	        static_cast<bool (*)(const Real&, const Real&, const Real&)>(&Eigen::internal::isApprox),
	        (py::arg("a"), "b", "eps"),
	        R"""(:return: ``bool``, True if ``a`` is approximately equal ``b`` with provided ``eps``, see also `here <https://stackoverflow.com/questions/15051367/how-to-compare-vectors-approximately-in-eigen>`__)""");
	py::def("isApproxOrLessThan",
	        static_cast<bool (*)(const Real&, const Real&, const Real&)>(&Eigen::internal::isApproxOrLessThan),
	        (py::arg("a"), "b", "eps"),
	        R"""(:return: ``bool``, True if ``a`` is approximately less than or equal ``b`` with provided ``eps``, see also `here <https://stackoverflow.com/questions/15051367/how-to-compare-vectors-approximately-in-eigen>`__)""");
#else
	// older eigen 3.2 didn't use `const Real&` but was copying third argument by value `Real`
	py::def("isMuchSmallerThan",
	        static_cast<bool (*)(const Real&, const Real&, Real)>(&Eigen::internal::isMuchSmallerThan<Real, Real>),
	        (py::arg("a"), "b", "eps"),
	        R"""(:return: ``bool``, True if ``a`` is less than ``b`` with provided ``eps``, see also `here <https://stackoverflow.com/questions/15051367/how-to-compare-vectors-approximately-in-eigen>`__)""");
	py::def("isApprox",
	        static_cast<bool (*)(const Real&, const Real&, Real)>(&Eigen::internal::isApprox<Real>),
	        (py::arg("a"), "b", "eps"),
	        R"""(:return: ``bool``, True if ``a`` is approximately equal ``b`` with provided ``eps``, see also `here <https://stackoverflow.com/questions/15051367/how-to-compare-vectors-approximately-in-eigen>`__)""");
	py::def("isApproxOrLessThan",
	        static_cast<bool (*)(const Real&, const Real&, Real)>(&Eigen::internal::isApproxOrLessThan<Real>),
	        (py::arg("a"), "b", "eps"),
	        R"""(:return: ``bool``, True if ``a`` is approximately less than or equal ``b`` with provided ``eps``, see also `here <https://stackoverflow.com/questions/15051367/how-to-compare-vectors-approximately-in-eigen>`__)""");
#endif

	py::def("toLongDouble",
	        static_cast<long double (*)(const Real&)>(&Eigen::internal::cast<Real, long double>),
	        (py::arg("x")),
	        R"""(:return: ``float`` converts ``Real`` type to ``long double`` and returns a native python ``float``.)""");
	py::def("toDouble",
	        static_cast<double (*)(const Real&)>(&Eigen::internal::cast<Real, double>),
	        (py::arg("x")),
	        R"""(:return: ``float`` converts ``Real`` type to ``double`` and returns a native python ``float``.)""");
	py::def("toLong",
	        static_cast<long (*)(const Real&)>(&Eigen::internal::cast<Real, long>),
	        (py::arg("x")),
	        R"""(:return: ``int`` converts ``Real`` type to ``long int`` and returns a native python ``int``.)""");
	py::def("toInt",
	        static_cast<int (*)(const Real&)>(&Eigen::internal::cast<Real, int>),
	        (py::arg("x")),
	        R"""(:return: ``int`` converts ``Real`` type to ``int`` and returns a native python ``int``.)""");

	expose_storage_ordering();

	// Very important line: Verfifies that Real type satisfies all the requirements of RealTypeConcept
	BOOST_CONCEPT_ASSERT((boost::math::concepts::RealTypeConcept<Real>));

	// check overload (and namespace) resolution for all math functions. As a side effect they are exported to python, and can be unit-tested.
#define YADE_PYEXPORT_MATH_1(func, DESC) py::def(#func, static_cast<Real (*)(const Real&)>(&::yade::math::func), (py::arg("x")), DESC);
#define YADE_PYEXPORT_MATH_1_MREF(func, DESC) py::def(#func, static_cast<Real (*)(Real)>(&::yade::math::func), (py::arg("x")), DESC);
#define YADE_PYEXPORT_MATH_1_COMPLEX(func, DESC) py::def(#func, static_cast<Complex (*)(const Complex&)>(&::yade::math::func), (py::arg("x")), DESC);
#define YADE_PYEXPORT_MATH_1_COMPLEX_TO_REAL(func, DESC) py::def(#func, static_cast<Real (*)(const Complex&)>(&::yade::math::func), (py::arg("x")), DESC);
#define YADE_PYEXPORT_MATH_1_INT(func, DESC) py::def(#func, static_cast<int (*)(const Real&)>(&::yade::math::func), (py::arg("x")), DESC);
	// FIXED: maybe registering the complex versions first will solve the problem that they are used when Real ones should be used - Yes.
	// Complex versions are registered first
	YADE_PYEXPORT_MATH_1_COMPLEX(
	        sin,
	        R"""(:return: ``Complex`` the sine of the ``Complex`` argument in radians. Depending on compilation options wraps ``::boost::multiprecision::sin(…)`` or `std::sin(…) <https://en.cppreference.com/w/cpp/numeric/complex/sin>`__ function.)""")
	YADE_PYEXPORT_MATH_1_COMPLEX(
	        sinh,
	        R"""(:return: ``Complex`` the hyperbolic sine of the ``Complex`` argument in radians. Depending on compilation options wraps ``::boost::multiprecision::sinh(…)`` or `std::sinh(…) <https://en.cppreference.com/w/cpp/numeric/complex/sinh>`__ function.)""")
	YADE_PYEXPORT_MATH_1_COMPLEX(
	        cos,
	        R"""(:return: ``Complex`` the cosine of the ``Complex`` argument in radians. Depending on compilation options wraps ``::boost::multiprecision::cos(…)`` or `std::cos(…) <https://en.cppreference.com/w/cpp/numeric/complex/cos>`__ function.)""")
	YADE_PYEXPORT_MATH_1_COMPLEX(
	        cosh,
	        R"""(:return: ``Complex`` the hyperbolic cosine of the ``Complex`` argument in radians. Depending on compilation options wraps ``::boost::multiprecision::cosh(…)`` or `std::cosh(…) <https://en.cppreference.com/w/cpp/numeric/complex/cosh>`__ function.)""")
	YADE_PYEXPORT_MATH_1_COMPLEX(
	        tan,
	        R"""(:return: ``Complex`` the tangent of the ``Complex`` argument in radians. Depending on compilation options wraps ``::boost::multiprecision::tan(…)`` or `std::tan(…) <https://en.cppreference.com/w/cpp/numeric/complex/tan>`__ function.)""")
	YADE_PYEXPORT_MATH_1_COMPLEX(
	        tanh,
	        R"""(:return: ``Complex`` the hyperbolic tangent of the ``Complex`` argument in radians. Depending on compilation options wraps ``::boost::multiprecision::tanh(…)`` or `std::tanh(…) <https://en.cppreference.com/w/cpp/numeric/complex/tanh>`__ function.)""")

	YADE_PYEXPORT_MATH_1_COMPLEX_TO_REAL(
	        abs,
	        R"""(:return: the ``Real`` absolute value of the ``Complex`` argument. Depending on compilation options wraps ``::boost::multiprecision::abs(…)`` or `std::abs(…) <https://en.cppreference.com/w/cpp/numeric/complex/abs>`__ function.)""")
	YADE_PYEXPORT_MATH_1_COMPLEX(
	        conj,
	        R"""(:return: the complex conjugation a ``Complex`` argument. Depending on compilation options wraps ``::boost::multiprecision::conj(…)`` or `std::conj(…) <https://en.cppreference.com/w/cpp/numeric/complex/conj>`__ function.)""")
	YADE_PYEXPORT_MATH_1_COMPLEX_TO_REAL(
	        real,
	        R"""(:return: the real part of a ``Complex`` argument. Depending on compilation options wraps ``::boost::multiprecision::real(…)`` or `std::real(…) <https://en.cppreference.com/w/cpp/numeric/complex/real2>`__ function.)""")
	YADE_PYEXPORT_MATH_1_COMPLEX_TO_REAL(
	        imag,
	        R"""(:return: the imag part of a ``Complex`` argument. Depending on compilation options wraps ``::boost::multiprecision::imag(…)`` or `std::imag(…) <https://en.cppreference.com/w/cpp/numeric/complex/imag2>`__ function.)""")

	YADE_PYEXPORT_MATH_1_COMPLEX(
	        exp,
	        R"""(:return: the base `e` exponential of a ``Complex`` argument. Depending on compilation options wraps ``::boost::multiprecision::exp(…)`` or `std::exp(…) <https://en.cppreference.com/w/cpp/numeric/complex/exp>`__ function.)""")
	YADE_PYEXPORT_MATH_1_COMPLEX(
	        log,
	        R"""(:return: the ``Complex`` natural (base `e`) logarithm of a complex value z with a branch cut along the negative real axis. Depending on compilation options wraps ``::boost::multiprecision::log(…)`` or `std::log(…) <https://en.cppreference.com/w/cpp/numeric/complex/log>`__ function.)""")

	// Real versions are registered afterwards
	YADE_PYEXPORT_MATH_1(
	        sin,
	        R"""(:return: ``Real`` the sine of the ``Real`` argument in radians. Depending on compilation options wraps ``::boost::multiprecision::sin(…)`` or `std::sin(…) <https://en.cppreference.com/w/cpp/numeric/math/sin>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        sinh,
	        R"""(:return: ``Real`` the hyperbolic sine of the ``Real`` argument in radians. Depending on compilation options wraps ``::boost::multiprecision::sinh(…)`` or `std::sinh(…) <https://en.cppreference.com/w/cpp/numeric/math/sinh>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        cos,
	        R"""(:return: ``Real`` the cosine of the ``Real`` argument in radians. Depending on compilation options wraps ``::boost::multiprecision::cos(…)`` or `std::cos(…) <https://en.cppreference.com/w/cpp/numeric/math/cos>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        cosh,
	        R"""(:return: ``Real`` the hyperbolic cosine of the ``Real`` argument in radians. Depending on compilation options wraps ``::boost::multiprecision::cosh(…)`` or `std::cosh(…) <https://en.cppreference.com/w/cpp/numeric/math/cosh>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        tan,
	        R"""(:return: ``Real`` the tangent of the ``Real`` argument in radians. Depending on compilation options wraps ``::boost::multiprecision::tan(…)`` or `std::tan(…) <https://en.cppreference.com/w/cpp/numeric/math/tan>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        tanh,
	        R"""(:return: ``Real`` the hyperbolic tangent of the ``Real`` argument in radians. Depending on compilation options wraps ``::boost::multiprecision::tanh(…)`` or `std::tanh(…) <https://en.cppreference.com/w/cpp/numeric/math/tanh>`__ function.)""")

#if (YADE_REAL_BIT > 64)
	YADE_PYEXPORT_MATH_1(
	        abs,
	        R"""(:return: the ``Real`` absolute value of the ``Real`` argument. Depending on compilation options wraps ``::boost::multiprecision::abs(…)`` or `std::abs(…) <https://en.cppreference.com/w/cpp/numeric/math/abs>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        fabs,
	        R"""(:return: the ``Real`` absolute value of the argument. Depending on compilation options wraps ``::boost::multiprecision::abs(…)`` or `std::abs(…) <https://en.cppreference.com/w/cpp/numeric/math/fabs>`__ function.)""")
#else
	YADE_PYEXPORT_MATH_1_MREF(
	        fabs,
	        R"""(:return: the ``Real`` absolute value of the argument. Depending on compilation options wraps ``::boost::multiprecision::abs(…)`` or `std::abs(…) <https://en.cppreference.com/w/cpp/numeric/math/fabs>`__ function.)""")
	YADE_PYEXPORT_MATH_1_MREF(
	        abs,
	        R"""(:return: the ``Real`` absolute value of the ``Real`` argument. Depending on compilation options wraps ``::boost::multiprecision::abs(…)`` or `std::abs(…) <https://en.cppreference.com/w/cpp/numeric/math/abs>`__ function.)""")
#endif
	YADE_PYEXPORT_MATH_1(
	        acos,
	        R"""(:return: ``Real`` the arcus cosine of the argument. Depending on compilation options wraps ``::boost::multiprecision::acos(…)`` or `std::acos(…) <https://en.cppreference.com/w/cpp/numeric/math/acos>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        acosh,
	        R"""(:return: ``Real`` the hyperbolic arcus cosine of the argument. Depending on compilation options wraps ``::boost::multiprecision::acosh(…)`` or `std::acosh(…) <https://en.cppreference.com/w/cpp/numeric/math/acosh>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        asin,
	        R"""(:return: ``Real`` the arcus sine of the argument. Depending on compilation options wraps ``::boost::multiprecision::asin(…)`` or `std::asin(…) <https://en.cppreference.com/w/cpp/numeric/math/asin>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        asinh,
	        R"""(:return: ``Real`` the hyperbolic arcus sine of the argument. Depending on compilation options wraps ``::boost::multiprecision::asinh(…)`` or `std::asinh(…) <https://en.cppreference.com/w/cpp/numeric/math/asinh>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        atan,
	        R"""(:return: ``Real`` the arcus tangent of the argument. Depending on compilation options wraps ``::boost::multiprecision::atan(…)`` or `std::atan(…) <https://en.cppreference.com/w/cpp/numeric/math/atan>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        atanh,
	        R"""(:return: ``Real`` the hyperbolic arcus tangent of the argument. Depending on compilation options wraps ``::boost::multiprecision::atanh(…)`` or `std::atanh(…) <https://en.cppreference.com/w/cpp/numeric/math/atanh>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        cbrt,
	        R"""(:return: ``Real`` cubic root of the argument. Depending on compilation options wraps ``::boost::multiprecision::cbrt(…)`` or `std::cbrt(…) <https://en.cppreference.com/w/cpp/numeric/math/cbrt>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        ceil,
	        R"""(:return: ``Real`` Computes the smallest integer value not less than arg. Depending on compilation options wraps ``::boost::multiprecision::ceil(…)`` or `std::ceil(…) <https://en.cppreference.com/w/cpp/numeric/math/ceil>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        erf,
	        R"""(:return: ``Real`` Computes the `error function <https://en.wikipedia.org/wiki/Error_function>`__ of argument. Depending on compilation options wraps ``::boost::multiprecision::erf(…)`` or `std::erf(…) <https://en.cppreference.com/w/cpp/numeric/math/erf>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        erfc,
	        R"""(:return: ``Real`` Computes the `complementary error function <https://en.wikipedia.org/wiki/Error_function#Complementary_error_function>`__ of argument, that is ``1.0-erf(arg)``. Depending on compilation options wraps ``::boost::multiprecision::erfc(…)`` or `std::erfc(…) <https://en.cppreference.com/w/cpp/numeric/math/erfc>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        exp,
	        R"""(:return: the base `e` exponential of a ``Real`` argument. Depending on compilation options wraps ``::boost::multiprecision::exp(…)`` or `std::exp(…) <https://en.cppreference.com/w/cpp/numeric/math/exp>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        exp2,
	        R"""(:return: the base `2` exponential of a ``Real`` argument. Depending on compilation options wraps ``::boost::multiprecision::exp2(…)`` or `std::exp2(…) <https://en.cppreference.com/w/cpp/numeric/math/exp2>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        expm1,
	        R"""(:return: the base `e` exponential of a ``Real`` argument minus ``1.0``. Depending on compilation options wraps ``::boost::multiprecision::expm1(…)`` or `std::expm1(…) <https://en.cppreference.com/w/cpp/numeric/math/expm1>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        floor,
	        R"""(:return: ``Real`` Computes the largest integer value not greater than arg. Depending on compilation options wraps ``::boost::multiprecision::floor(…)`` or `std::floor(…) <https://en.cppreference.com/w/cpp/numeric/math/floor>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        ilogb,
	        R"""(:return: ``Real`` extracts the value of the unbiased exponent from the floating-point argument arg, and returns it as a signed integer value. Depending on compilation options wraps ``::boost::multiprecision::ilogb(…)`` or `std::ilogb(…) <https://en.cppreference.com/w/cpp/numeric/math/ilogb>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        lgamma,
	        R"""(:return: ``Real`` Computes the natural logarithm of the absolute value of the `gamma function <https://en.wikipedia.org/wiki/Gamma_function>`__ of arg. Depending on compilation options wraps ``::boost::multiprecision::lgamma(…)`` or `std::lgamma(…) <https://en.cppreference.com/w/cpp/numeric/math/lgamma>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        log,
	        R"""(:return: the ``Real`` natural (base `e`) logarithm of a real value. Depending on compilation options wraps ``::boost::multiprecision::log(…)`` or `std::log(…) <https://en.cppreference.com/w/cpp/numeric/math/log>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        log10,
	        R"""(:return: the ``Real`` decimal (base ``10``) logarithm of a real value. Depending on compilation options wraps ``::boost::multiprecision::log10(…)`` or `std::log10(…) <https://en.cppreference.com/w/cpp/numeric/math/log10>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        log1p,
	        R"""(:return: the ``Real`` natural (base `e`) logarithm of ``1+argument``. Depending on compilation options wraps ``::boost::multiprecision::log1p(…)`` or `std::log1p(…) <https://en.cppreference.com/w/cpp/numeric/math/log1p>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        log2,
	        R"""(:return: the ``Real`` binary (base ``2``) logarithm of a real value. Depending on compilation options wraps ``::boost::multiprecision::log2(…)`` or `std::log2(…) <https://en.cppreference.com/w/cpp/numeric/math/log2>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        logb,
	        R"""(:return: Extracts the value of the unbiased radix-independent exponent from the floating-point argument arg, and returns it as a floating-point value. Depending on compilation options wraps ``::boost::multiprecision::logb(…)`` or `std::logb(…) <https://en.cppreference.com/w/cpp/numeric/math/logb>`__ function.)""")
	//YADE_PYEXPORT_MATH_1(riemann_zeta) // since C++17
	YADE_PYEXPORT_MATH_1(
	        rint,
	        R"""(:return: Rounds the floating-point argument arg to an integer value (in floating-point format), using the `current rounding mode <https://en.cppreference.com/w/cpp/numeric/fenv/FE_round>`__. Depending on compilation options wraps ``::boost::multiprecision::rint(…)`` or `std::rint(…) <https://en.cppreference.com/w/cpp/numeric/math/rint>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        round,
	        R"""(:return: ``Real`` the nearest integer value to arg (in floating-point format), rounding halfway cases away from zero, regardless of the current rounding mode.. Depending on compilation options wraps ``::boost::multiprecision::round(…)`` or `std::round(…) <https://en.cppreference.com/w/cpp/numeric/math/round>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        sqrt,
	        R"""(:return: ``Real`` square root of the argument. Depending on compilation options wraps ``::boost::multiprecision::sqrt(…)`` or `std::sqrt(…) <https://en.cppreference.com/w/cpp/numeric/math/sqrt>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        tgamma,
	        R"""(:return: ``Real`` Computes the `gamma function <https://en.wikipedia.org/wiki/Gamma_function>`__ of arg. Depending on compilation options wraps ``::boost::multiprecision::tgamma(…)`` or `std::tgamma(…) <https://en.cppreference.com/w/cpp/numeric/math/tgamma>`__ function.)""")
	YADE_PYEXPORT_MATH_1(
	        trunc,
	        R"""(:return: ``Real`` the nearest integer not greater in magnitude than arg. Depending on compilation options wraps ``::boost::multiprecision::trunc(…)`` or `std::trunc(…) <https://en.cppreference.com/w/cpp/numeric/math/trunc>`__ function.)""")

	YADE_PYEXPORT_MATH_1_INT(sgn, R"""(:return: ``int`` the sign of the argument: ``-1``, ``0`` or ``1``.)""")
	YADE_PYEXPORT_MATH_1_INT(sign, R"""(:return: ``int`` the sign of the argument: ``-1``, ``0`` or ``1``.)""")
#undef YADE_PYEXPORT_MATH_1
#undef YADE_PYEXPORT_MATH_1_MREF
#undef YADE_PYEXPORT_MATH_1_COMPLEX
#undef YADE_PYEXPORT_MATH_1_COMPLEX_TO_REAL
#undef YADE_PYEXPORT_MATH_1_INT

#define YADE_PYEXPORT_MATH_2(func, DESC) py::def(#func, static_cast<Real (*)(const Real&, const Real&)>(&::yade::math::func), (py::arg("x"), "y"), DESC);
#define YADE_PYEXPORT_MATH_2_CREF(func, DESC)                                                                                                                  \
	py::def(#func,                                                                                                                                         \
	        static_cast<const Real& (*)(const Real&, const Real&)>(&::yade::math::func),                                                                   \
	        (py::arg("x"), "y"),                                                                                                                           \
	        py::return_value_policy<py::copy_const_reference>(),                                                                                           \
	        DESC);
	YADE_PYEXPORT_MATH_2(
	        atan2,
	        R"""(:return: ``Real`` the arc tangent of y/x using the signs of the arguments ``x`` and ``y`` to determine the correct quadrant. Depending on compilation options wraps ``::boost::multiprecision::atan2(…)`` or `std::atan2(…) <https://en.cppreference.com/w/cpp/numeric/math/atan2>`__ function.)""")
	//YADE_PYEXPORT_MATH_2(beta) // since C++17
	//YADE_PYEXPORT_MATH_2(cyl_bessel_i) // since C++17
	//YADE_PYEXPORT_MATH_2(cyl_bessel_j) // since C++17
	//YADE_PYEXPORT_MATH_2(cyl_bessel_k) // since C++17
	YADE_PYEXPORT_MATH_2(
	        fmod,
	        R"""(:return: ``Real`` the floating-point remainder of the division operation ``x/y`` of the arguments ``x`` and ``y``. Depending on compilation options wraps ``::boost::multiprecision::fmod(…)`` or `std::fmod(…) <https://en.cppreference.com/w/cpp/numeric/math/fmod>`__ function.)""")
	YADE_PYEXPORT_MATH_2(
	        hypot,
	        R"""(:return: ``Real`` the square root of the sum of the squares of ``x`` and ``y``, without undue overflow or underflow at intermediate stages of the computation. Depending on compilation options wraps ``::boost::multiprecision::hypot(…)`` or `std::hypot(…) <https://en.cppreference.com/w/cpp/numeric/math/hypot>`__ function.)""")
	YADE_PYEXPORT_MATH_2_CREF(
	        max,
	        R"""(:return: ``Real`` larger of the two arguments. Depending on compilation options wraps ``::boost::multiprecision::max(…)`` or `std::max(…) <https://en.cppreference.com/w/cpp/numeric/math/max>`__ function.)""")
	YADE_PYEXPORT_MATH_2_CREF(
	        min,
	        R"""(:return: ``Real`` smaller of the two arguments. Depending on compilation options wraps ``::boost::multiprecision::min(…)`` or `std::min(…) <https://en.cppreference.com/w/cpp/numeric/math/min>`__ function.)""")
	YADE_PYEXPORT_MATH_2(
	        pow,
	        R"""(:return: ``Real`` the value of ``base`` raised to the power ``exp``. Depending on compilation options wraps ``::boost::multiprecision::pow(…)`` or `std::pow(…) <https://en.cppreference.com/w/cpp/numeric/math/pow>`__ function.)""")
	YADE_PYEXPORT_MATH_2(
	        remainder,
	        R"""(:return: ``Real`` the IEEE remainder of the floating point division operation ``x/y``. Depending on compilation options wraps ``::boost::multiprecision::remainder(…)`` or `std::remainder(…) <https://en.cppreference.com/w/cpp/numeric/math/remainder>`__ function.)""")
#undef YADE_PYEXPORT_MATH_2
#undef YADE_PYEXPORT_MATH_2_CREF

#define YADE_PYEXPORT_MATH_2_TYPE1(func, FirstType) py::def(#func, static_cast<Real (*)(FirstType, const Real&)>(&::yade::math::func), (py::arg("x"), "y"));
	//YADE_PYEXPORT_MATH_2_TYPE1(sph_bessel, unsigned) // since C++17
#undef YADE_PYEXPORT_MATH_2_TYPE1

#define YADE_PYEXPORT_MATH_2_TYPE2(func, SecondType, DESC)                                                                                                     \
	py::def(#func, static_cast<Real (*)(const Real&, SecondType)>(&::yade::math::func), (py::arg("x"), "y"), DESC);
	YADE_PYEXPORT_MATH_2_TYPE2(
	        ldexp,
	        int,
	        R"""(:return: Multiplies a floating point value ``x`` by the number 2 raised to the ``exp`` power. Depending on compilation options wraps ``::boost::multiprecision::ldexp(…)`` or `std::ldexp(…) <https://en.cppreference.com/w/cpp/numeric/math/ldexp>`__ function.)""")
#undef YADE_PYEXPORT_MATH_2_TYPE2

	std_pair_to_python_converter<Real, Real>();
	std_pair_to_python_converter<double, double>();
	std_pair_to_python_converter<Real, long>();
	std_pair_to_python_converter<Real, int>();
	py::def("frexp",
	        test_frexp,
	        (py::arg("x")),
	        R"""(:return: tuple of ``(Real,int)``, decomposes given floating point ``Real`` argument into a normalized fraction and an integral power of two. Depending on compilation options wraps ``::boost::multiprecision::frexp(…)`` or `std::frexp(…) <https://en.cppreference.com/w/cpp/numeric/math/frexp>`__ function.)""");
	py::def("modf",
	        test_modf,
	        (py::arg("x")),
	        R"""(:return: tuple of ``(Real,Real)``, decomposes given floating point ``Real`` into integral and fractional parts, each having the same type and sign as x. Depending on compilation options wraps ``::boost::multiprecision::modf(…)`` or `std::modf(…) <https://en.cppreference.com/w/cpp/numeric/math/modf>`__ function.)""");

#define YADE_PYEXPORT_MATH_3(func, DESC)                                                                                                                       \
	py::def(#func, static_cast<Real (*)(const Real&, const Real&, const Real&)>(&::yade::math::func), (py::arg("x"), "y", "z"), DESC);
	YADE_PYEXPORT_MATH_3(
	        fma,
	        R"""(:return: ``Real`` - computes ``(x*y) + z`` as if to infinite precision and rounded only once to fit the result type. Depending on compilation options wraps ``::boost::multiprecision::fma(…)`` or `std::fma(…) <https://en.cppreference.com/w/cpp/numeric/math/fma>`__ function.)""")
#undef YADE_PYEXPORT_MATH_3

	py::def("remquo",
	        test_remquo,
	        (py::arg("x"), "y"),
	        R"""(:return: tuple of ``(Real,long)``, the floating-point remainder of the division operation ``x/y`` as the std::remainder() function does. Additionally, the sign and at least the three of the last bits of ``x/y`` are returned, sufficient to determine the octant of the result within a period. Depending on compilation options wraps ``::boost::multiprecision::remquo(…)`` or `std::remquo(…) <https://en.cppreference.com/w/cpp/numeric/math/remquo>`__ function.)""");
	py::def("testArray", ::yade::testArray, R"""(This function tests call to ``std::vector::data(…)`` function in order to extract the array.)""");

#ifdef YADE_CGAL
	py::scope().attr("testCgalNumTraits") = true;
	// https://doc.cgal.org/latest/Algebraic_foundations/group__PkgAlgebraicFoundationsRef.html
	py::def("CGAL_Is_valid",
	        test_CGAL_Is_valid,
	        (py::arg("x")),
	        R"""(
CGAL's function ``Is_valid``, as described in `CGAL algebraic <https://doc.cgal.org/latest/Algebraic_foundations/index.html>`__
`foundations <https://doc.cgal.org/latest/Algebraic_foundations/group__PkgAlgebraicFoundationsRef.html>`__ :ysrc:`exposed<lib/high-precision/CgalNumTraits.hpp>`
to python for :ysrccommit:`testing<ff600a80018d21c03626c720cda08967b043c1c8/py/tests/testMath.py#L207>` of CGAL numerical traits.

:return: ``bool`` indicating if the ``Real`` argument is valid. Checks are performed against NaN and Inf.
)""");
	// AlgebraicStructureTraits
	py::def("CGAL_Square",
	        test_CGAL_Square,
	        (py::arg("x")),
	        R"""(
CGAL's function ``Square``, as described in `CGAL algebraic <https://doc.cgal.org/latest/Algebraic_foundations/index.html>`__
`foundations <https://doc.cgal.org/latest/Algebraic_foundations/group__PkgAlgebraicFoundationsRef.html>`__ :ysrc:`exposed<lib/high-precision/CgalNumTraits.hpp>`
to python for :ysrccommit:`testing<ff600a80018d21c03626c720cda08967b043c1c8/py/tests/testMath.py#L207>` of CGAL numerical traits.

:return: ``Real`` the argument squared.
)""");
	py::def("CGAL_Sqrt",
	        test_CGAL_Sqrt,
	        (py::arg("x")),
	        R"""(
CGAL's function ``Sqrt``, as described in `CGAL algebraic <https://doc.cgal.org/latest/Algebraic_foundations/index.html>`__
`foundations <https://doc.cgal.org/latest/Algebraic_foundations/group__PkgAlgebraicFoundationsRef.html>`__ :ysrc:`exposed<lib/high-precision/CgalNumTraits.hpp>`
to python for :ysrccommit:`testing<ff600a80018d21c03626c720cda08967b043c1c8/py/tests/testMath.py#L207>` of CGAL numerical traits.

:return: ``Real`` the square root of argument.
)""");
	py::def("CGAL_Kth_root",
	        test_CGAL_Kth_root,
	        (py::arg("x")),
	        R"""(
CGAL's function ``Kth_root``, as described in `CGAL algebraic <https://doc.cgal.org/latest/Algebraic_foundations/index.html>`__
`foundations <https://doc.cgal.org/latest/Algebraic_foundations/group__PkgAlgebraicFoundationsRef.html>`__ :ysrc:`exposed<lib/high-precision/CgalNumTraits.hpp>`
to python for :ysrccommit:`testing<ff600a80018d21c03626c720cda08967b043c1c8/py/tests/testMath.py#L207>` of CGAL numerical traits.

:return: ``Real`` the k-th root of argument.
)""");
	// RealEmbeddableTraits
	py::def("CGAL_To_interval",
	        test_CGAL_To_interval,
	        (py::arg("x")),
	        R"""(
CGAL's function ``To_interval``, as described in `CGAL algebraic <https://doc.cgal.org/latest/Algebraic_foundations/index.html>`__
`foundations <https://doc.cgal.org/latest/Algebraic_foundations/group__PkgAlgebraicFoundationsRef.html>`__ :ysrc:`exposed<lib/high-precision/CgalNumTraits.hpp>`
to python for :ysrccommit:`testing<ff600a80018d21c03626c720cda08967b043c1c8/py/tests/testMath.py#L207>` of CGAL numerical traits.

:return: ``(double,double)`` tuple inside which the high-precision ``Real`` argument resides.
)""");
	py::def("CGAL_Sgn",
	        test_CGAL_Sgn,
	        (py::arg("x")),
	        R"""(
CGAL's function ``Sgn``, as described in `CGAL algebraic <https://doc.cgal.org/latest/Algebraic_foundations/index.html>`__
`foundations <https://doc.cgal.org/latest/Algebraic_foundations/group__PkgAlgebraicFoundationsRef.html>`__ :ysrc:`exposed<lib/high-precision/CgalNumTraits.hpp>`
to python for :ysrccommit:`testing<ff600a80018d21c03626c720cda08967b043c1c8/py/tests/testMath.py#L207>` of CGAL numerical traits.

:return: sign of the argument, can be ``-1``, ``0`` or ``1``. Not very useful in python. In C++ it is useful to obtain a sign of an expression with exact accuracy, CGAL starts using MPFR internally for this when the approximate interval contains zero inside it.
)""");
	py::def("CGAL_Is_finite",
	        test_CGAL_Is_finite,
	        (py::arg("x")),
	        R"""(
CGAL's function ``Is_finite``, as described in `CGAL algebraic <https://doc.cgal.org/latest/Algebraic_foundations/index.html>`__
`foundations <https://doc.cgal.org/latest/Algebraic_foundations/group__PkgAlgebraicFoundationsRef.html>`__ :ysrc:`exposed<lib/high-precision/CgalNumTraits.hpp>`
to python for :ysrccommit:`testing<ff600a80018d21c03626c720cda08967b043c1c8/py/tests/testMath.py#L207>` of CGAL numerical traits.

:return: ``bool`` indicating if the ``Real`` argument is finite.
)""");
	py::def("CGAL_simpleTest",
	        ::yade::testSimpleCgalNumTraitsCalculation,
	        R"""(
Tests a simple CGAL calculation. Distance between plane and point, uses CGAL's sqrt and pow.

:return: 3.0
)""");
#else
	py::scope().attr("testCgalNumTraits") = false;
#endif

} catch (...) {
	LOG_FATAL("Importing this module caused an exception and this module is in an inconsistent state now.");
	PyErr_Print();
	PyErr_SetString(PyExc_SystemError, __FILE__); // raising anything other than SystemError is not possible
	boost::python::handle_exception();
	throw;
}

