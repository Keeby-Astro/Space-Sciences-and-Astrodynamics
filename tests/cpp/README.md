# C++ Test Scripts

This directory contains unit and validation tests for the C++ components of Space-Sciences-and-Astrodynamics.

## Structure
- All `.cpp` files in this directory are built as test executables by CMake.
- Tests are linked against the main library and, if available, Catch2 for modern C++ testing.

## Running Tests
1. Configure and build the project with CMake:
   ```powershell
   cmake -S . -B build
   cmake --build build
   ```
2. Run all tests (if Catch2 is found):
   ```powershell
   ctest --test-dir build
   ```
   Or run individual test executables directly from the build directory.

## Adding Tests
- Add new `.cpp` files for each test case or suite.
- Use [Catch2](https://github.com/catchorg/Catch2) for expressive, header-only test cases:
   ```cpp
   #define CATCH_CONFIG_MAIN
   #include <catch2/catch.hpp>
   #include "sv_from_coe.h"

   TEST_CASE("sv_from_coe computes correct state vector", "[orbit]") {
       // ... test code ...
   }
   ```
- Tests can also be written using other frameworks or plain asserts if needed.

---

For more information, see the main project README and CMakeLists.txt.
