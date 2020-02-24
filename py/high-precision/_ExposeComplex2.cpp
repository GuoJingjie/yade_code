#ifdef _HIGH_PRECISION_SUPPORT
// wall clock time: 0:24.22 → split into two files
#include <lib/high-precision/Real.hpp>
#include <lib/high-precision/ToFromPythonConverter.hpp>
using namespace ::yade::MathEigenTypes;
// define this for compatibility with minieigen.
#define _COMPLEX_SUPPORT
#include <minieigen/expose-complex.cpp>

#endif

