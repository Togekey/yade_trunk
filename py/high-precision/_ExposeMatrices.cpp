#include <lib/high-precision/Real.hpp>
#include <lib/high-precision/ToFromPythonConverter.hpp>
using namespace ::yade::MathEigenTypes;
#ifdef EXPERIMENTS_ONLY_LOCAL_MINIEIGEN
#include <minieigen-local/expose-matrices.cpp>
#else
#include <minieigen/expose-matrices.cpp>
#endif

