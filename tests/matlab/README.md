# MATLAB Test Scripts

This directory contains unit and validation tests for the MATLAB components of Space-Sciences-and-Astrodynamics.

## Structure
- Place all MATLAB `.m` test scripts and functions here.
- Organize tests using MATLAB's [unit testing framework](https://www.mathworks.com/help/matlab/matlab-unit-test-framework.html) when possible.
- Tests should cover all major algorithms, edge cases, and numerical methods implemented in the `curtis_scripts/matlab` namespace.

## Running Tests
1. Open MATLAB and set the working directory to this folder:
   ```matlab
   cd('.../Space-Sciences-and-Astrodynamics/tests/matlab')
   ```
2. Run all tests using:
   ```matlab
   runtests
   ```
   Or run individual test files as needed.

## Adding Tests
- Create new `.m` files for each test case or suite.
- Use the `matlab.unittest.TestCase` class for structured, automated tests:
   ```matlab
   classdef test_sv_from_coe < matlab.unittest.TestCase
       methods(Test)
           function testBasicCase(tc)
               % ... test code ...
           end
       end
   end
   ```

---

For more information, see the main project README and MATLAB documentation.
