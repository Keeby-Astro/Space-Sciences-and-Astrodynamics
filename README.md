# Space-Sciences-and-Astrodynamics

This repository contains MATLAB and Python implementations of several orbital mechanics and astrodynamics functions. The toolkit offers functionality for simulating orbits, computing transformations, performing numerical integrations, and more.

## Table of Contents
- [Overview](#overview)
- [Requirements](#requirements)
- [Features](#features)
- [Files and Modules](#files-and-modules)
- [Usage](#usage)
- [Installation](#installation)

---

## Overview
The **Space-Sciences-and-Astrodynamics** Repository provides a set of MATLAB and Python scripts based on Orbital Mechanics for Engineering Students by Howard D. Curtis for modeling and simulating space-based mechanics, including coordinate transformations, orbital element calculations, interplanetary trajectory solutions, and more. The toolkit is designed for educational and research purposes and is a resource for those interested in understanding or working within astrodynamics.

## Requirements
- **MATLAB** (R2021a or later recommended) OR **Python** 3.8+
- MATLAB’s `ode45` and other built-in functions
- Python packages (only if using Python versions of the scripts):
  - `NumPy`
  - `SciPy`
  - `Matplotlib`

## Features
- **Coordinate Transformation**: Convert between orbital elements and state vectors.
- **Orbit Simulations**: Simulate different types of orbits, including two-body and three-body systems.
- **Trajectory Calculation**: Calculate interplanetary and lunar trajectories.
- **Root Finding & Numerical Methods**: Implement bisection, Runge-Kutta, and Heun's methods for solving differential equations.
- **Orbital Elements & Ephemerides**: Compute positions of celestial bodies (e.g., Sun, Moon) for a given Julian date.

## Files and Modules

### Course Code Modules
| File                      | Description |
|---------------------------|-------------|
| `Coordinate Transformation.py` | Converts between orbital elements and state vectors for orbital dynamics, calculating angular momentum, eccentricity, inclination, and other essential values. |
| `Family_of_orbits.py`         | Simulates and visualizes a family of orbits with varying semi-major axis and eccentricity. Also includes an animation function for orbit visualization. |
| `Orbital Elements.py`         | Calculates orbital parameters from initial state vectors, such as vector magnitude, radial velocity, and angular momentum. |
| `Orbit_Simulator.py`          | Simulates orbital motion with transformations between orbital elements and Cartesian coordinates, providing a 3D orbit visualization. |

### MATLAB and Python Code Modules
#### `atan2d_0_360.m` | `atan2d_0_360.py`

- The `atan2d_0_360` function is a specialized variant of the traditional `atan2` function, designed for scenarios where:
  - Angles need to be normalized to the range [0, 360] instead of the standard [-180, 180].
  - Consistent angle representation is required for astrodynamics, robotics, and other applications.
  - Improved readability and usability of angular calculations are beneficial in workflows.

- The test scripts included in both MATLAB and Python versions validate the `atan2d_0_360` function across a variety of input cases. These scripts:
  - Test the function with edge cases such as zero and axis-aligned vectors.
  - Evaluate performance on quadrants to ensure correct angle normalization.
  - Demonstrate the function's application in common scenarios.

![atan2d_0_360_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/atan2d_0_360_plot.png)

## Usage

### MATLAB
Clone the repository and navigate to the MATLAB files directory.

### Python
Clone the repository and navigate to the Python files directory.

## Installation
```bash
git clone https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics.git
```

## Future Updates
C++ and Fortran translations are currently in progress. This READme file will constantly be updated as needed.
