# Python Test Scripts

This directory contains unit and validation tests for the Python components of Space-Sciences-and-Astrodynamics.

## Structure
- `unit/` — Unit tests for individual functions and modules.
- `validation/` — Validation tests comparing outputs to reference data or published results.

## Running Tests
To run all tests, navigate to the root of the repository and use:

```powershell
python -m unittest discover tests/python
```
Or, if using `pytest`:
```powershell
pytest tests/python
```
