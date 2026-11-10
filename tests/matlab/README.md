# MATLAB Test Scripts

This directory contains runnable examples and validation tests for the MATLAB components of Space-Sciences-and-Astrodynamics. Some files print results or draw plots; a few also use assertions.

## Structure
- Place all MATLAB `.m` test scripts and functions here.
- Organize tests using MATLAB's [unit testing framework](https://www.mathworks.com/help/matlab/matlab-unit-test-framework.html) when possible.
- Tests should cover all major algorithms, edge cases, and numerical methods implemented in the `curtis_scripts/matlab` folder.

## Running Tests
Open a test file in the MATLAB Editor and click **Run**. If MATLAB asks, change to the test's folder or add that folder to the path. You can also make `tests/matlab` the current folder and run a test by name:

```matlab
atan2d_0_360_test
stumpS_test
```

Each test adds `curtis_scripts/matlab` to the search path for the current MATLAB session. The path is calculated from the test file's location, so moving the repository or using another operating system does not require editing the path. Keep the repository's folder structure intact and run the whole file so this setup executes.

From the repository root, you can run a script test with:

```matlab
run(fullfile('tests', 'matlab', 'atan2d_0_360_test.m'))
```

To call the function-based examples from the repository root, first make their folder visible:

```matlab
addpath(fullfile(pwd, 'tests', 'matlab'))
stumpS_test
```

Run these examples individually. The existing files are ordinary scripts and callable functions, rather than suites written for MATLAB's `runtests` framework. Completing without an error does not verify numerical accuracy for examples that only print results or plot.

If a function is still unrecognized, use `which function_name -all` to check whether MATLAB can find it. Repository helpers must be present in `curtis_scripts/matlab`; external dependencies must be installed separately. For example, `quat_rotate` requires the Aerospace Toolbox functions `quatmultiply` and `quatinv`. These path statements do not install toolboxes or save permanent changes to your MATLAB path.

## Adding Tests
- Create new `.m` files for each test case or suite.
- For another directly runnable example in this folder, put this statement after the introductory comments (inside the main function, if there is one), before calling repository helpers:
   ```matlab
   addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', ...
       'curtis_scripts', 'matlab'));
   ```
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
