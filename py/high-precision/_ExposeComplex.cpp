#include <lib/high-precision/Real.hpp>
#include <lib/high-precision/ToFromPythonConverter.hpp>
using namespace ::yade::MathEigenTypes;
// define this for compatibility with minieigen.
#define _COMPLEX_SUPPORT
#ifdef EXPERIMENTS_ONLY_LOCAL_MINIEIGEN
#include <minieigen-local/expose-complex.cpp>
#else
#include <minieigen/expose-complex.cpp>
#endif

