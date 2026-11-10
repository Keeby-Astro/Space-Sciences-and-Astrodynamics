Script Reference
================

The repository contains standalone educational scripts. MATLAB, C++ and Python
implementations retain their existing language-specific organization and interfaces.

Python Functions
----------------

Add ``curtis_scripts/python`` to the Python import path, or install the flat
modules using ``python -m pip install .`` from the repository root. Import the
root finder with ``from bisect_root import bisect``.

.. automodule:: coe_from_sv
   :members:

.. automodule:: sv_from_coe
   :members:

.. automodule:: rv_from_r0v0
   :members:

.. automodule:: lambert
   :members:

.. automodule:: rkf45
   :members:

C++ Scripts
-----------

The programs in ``curtis_scripts/cpp/sources`` compile individually using the
headers in ``curtis_scripts/cpp/headers``. CMake builds each program as a separate
executable. Python bindings are currently a placeholder.

MATLAB Scripts
--------------

Add ``curtis_scripts/matlab`` to the MATLAB path. Functions that declare global
parameters require those values to be initialized by the caller. Several examples
use MATLAB toolbox functions; consult the audit report for dependency details.

Course Examples
---------------

The four demonstrations in ``course_scripts`` retain their course-specific
orbital-element ordering. ``functions.py`` supplies the coordinate conversions
used by the orbit simulator.

Algorithm Mapping
-----------------

The repository-root ``Algorithm_Number_Reference.md`` maps the textbook algorithms
to their source files. The audit report records numerical corrections and validation.
