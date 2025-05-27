# C++ Python Bindings (pybind11)

This directory contains pybind11-based bindings for exposing C++ orbital mechanics and astrodynamics routines to Python. These bindings enable high-performance C++ code to be used directly from Python, supporting both research and educational workflows.

## Structure
- Place all binding source files (e.g., `bindings.cpp`, `module.cpp`) here.
- Bindings should wrap C++ functions, classes, and data structures from the main library in `curtis_scripts/cpp/Sources`.
- The resulting Python extension module will be built as part of the monorepo using `scikit-build-core` and CMake.

## Getting Started
- See the [pybind11 documentation](https://pybind11.readthedocs.io/) for binding patterns and best practices.
- Example CMake integration is provided in the root `CMakeLists.txt`.

## Example
```cpp
#include <pybind11/pybind11.h>
#include "sv_from_coe.h"

PYBIND11_MODULE(space_mech_cpp, m) {
    m.def("sv_from_coe", &sv_from_coe, "State vector from classical orbital elements");
}
```