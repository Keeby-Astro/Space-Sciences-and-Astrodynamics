# Orbital Mechanics Simulation Toolkit

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
| `twobody3d.m`                 | Solves the two-body gravitational problem in 3D using Runge-Kutta-Fehlberg 4(5), calculating trajectories for both masses and their center of mass. |
| `threebody.m`                 | Simulates the planar three-body problem using MATLAB’s `ode45` solver, showing trajectories relative to the inertial frame and center of mass. |

### MATLAB Code Modules
| File                     | Description |
|--------------------------|-------------|
| `atan2d_0_360.m`         | Computes the arctangent of \( y/x \) in degrees, constrained to [0, 360]. |
| `atmosphere.m`           | Computes atmospheric density up to 1000 km using exponential interpolation for various atmospheric layers. |
| `bisect.m`               | Implements the bisection method to find roots within specified intervals. |
| `dcm_to_euler.m`, `dcm_to_ypr.m`, `dcm_from_q.m` | Direction cosine matrix (DCM) transformations to Euler angles and quaternion formats. |
| `f_and_g.m`, `fDot_and_gDot.m`, `f_and_g_ta.m` | Calculates Lagrange coefficients for orbital propagation over time or with a change in true anomaly. |
| `kepler_E.m`, `kepler_H.m`, `kepler_U.m` | Solves Kepler's equations for elliptical, hyperbolic, and universal anomaly cases. |
| `lambert.m`              | Solves Lambert's problem for orbital transfer calculations between two points in space. |
| `los.m`                  | Determines whether the Earth obstructs the line of sight between a satellite and the Sun. |
| `LST.m`                  | Calculates local sidereal time based on location and universal time. |
| `lunar_position.m`       | Computes the Moon's geocentric equatorial position for a given Julian date. |
| `month_planet_names.m`   | Maps month and planet IDs to their respective names. |
| `planet_elements_and_sv.m` | Calculates a planet's orbital elements and state vector for a given date and time. |
| `q_from_dcm.m`           | Converts a DCM into a quaternion representation for coordinate transformations. |
| `ra_and_dec_from_r.m`    | Converts a position vector to right ascension and declination angles. |
| `rk1_4.m`, `rkf45.m`     | Implements various Runge-Kutta methods (RK1-RK4, RKF4(5)) for differential equation integration. |
| `rv_from_observe.m`      | Calculates position and velocity from radar observations. |
| `rv_from_r0v0.m`, `rv_from_r0v0_ta.m` | Propagates a state vector forward in time or through a change in true anomaly. |
| `rva_relative.m`         | Calculates relative position, velocity, and acceleration of two objects in an LVLH frame. |
| `simpsons_lunar_ephemeris.m` | Computes the Moon's position and velocity vectors based on Simpson's ephemeris model. |
| `solar_position.m`       | Calculates the position vector of the Sun relative to Earth for a given Julian date. |
| `stumpC.m`, `stumpS.m`   | Evaluates Stumpff functions, used in solving Kepler's equations. |

### Python Code Modules
Python versions of all above modules are also available, with the same functionality, optimized for `NumPy` and `SciPy`.

## Usage

### MATLAB
Clone the repository and navigate to the MATLAB files directory.

### Python
Clone the repository and navigate to the Python files directory.

## Installation
```bash
git clone https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics.git
