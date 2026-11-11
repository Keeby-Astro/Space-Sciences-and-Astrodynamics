# Python Test Scripts

This directory contains directly runnable examples and validation scripts for the Python components of Space-Sciences-and-Astrodynamics. Most print results or draw plots; `stumpC_test.py` and `stumpS_test.py` also use assertions.

## Running a Script

Install the external dependencies in the Python environment you use to run the scripts:

```sh
python -m pip install numpy scipy matplotlib
```

From the repository root, run an individual script:

```sh
python tests/python/atan2d_0_360_test.py
python tests/python/stumpS_test.py
```

If your current folder is `tests/python`, use:

```sh
python atan2d_0_360_test.py
```

You can also use your editor's command to run the whole Python file with the same interpreter. From any other working folder, supply the script's absolute path, quoting paths that contain spaces.

Each test locates `curtis_scripts/python` relative to its own `__file__` and adds that folder to the current Python process's module search path. Keep the repository's folder structure intact. You do not need to install this repository as a package, set `PYTHONPATH`, or edit machine-specific paths. External dependencies still need to be installed separately.

The source folder is appended to the search path to preserve standard-library imports. `bisect_test.py` explicitly loads the repository's `bisect.py`, whose root-finding function shares a module name with Python's standard-library `bisect`. This leaves the standard-library module available to Python and external dependencies.

These files use ordinary functions and `if __name__ == "__main__"` blocks. Run them individually; `pytest` and `unittest` discovery do not execute those main blocks. Completing without an error does not verify numerical accuracy for examples that only print results or plot. Plotting examples open figure windows when run with a graphical backend.

## Adding Scripts

For another directly runnable test in this folder, place this block before importing repository helpers:

```python
import sys
from pathlib import Path

_SCRIPTS_DIR = str(Path(__file__).resolve().parents[2] / "curtis_scripts" / "python")
if _SCRIPTS_DIR not in sys.path:
    sys.path.append(_SCRIPTS_DIR)
```

## Troubleshooting Imports

For `ModuleNotFoundError`, check whether the missing module is a repository helper or an external dependency. Repository helpers must be present in `curtis_scripts/python`; install external dependencies with the same interpreter that runs the script.

If an import finds an unexpected module, inspect its location after the import, for example:

```python
print(atan2d_0_360.__file__)
```

An installed module with the same name can take precedence over the appended source folder. Run in an environment without that conflicting installation. Restart an interactive session after switching between repository copies, since Python caches imported modules.
