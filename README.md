# Space-Sciences-and-Astrodynamics

This repository contains MATLAB, Python, and C++ implementations of orbital mechanics and astrodynamics tools. It provides a comprehensive set of functions for simulating orbital motion, performing numerical methods, analyzing spacecraft trajectories, and modeling celestial mechanics. The repository is based on foundational principles and extends concepts found in *Orbital Mechanics for Engineering Students* by Howard D. Curtis.

## Table of Contents
- [Overview](#overview)
- [Requirements](#requirements)
- [Features](#features)
- [Course Code Modules](#course-code-modules)
- [MATLAB Python and C++ Code Modules](#matlab-python-and-c-code-modules)
- [Usage](#usage)
- [Installation](#installation)
- [Future Updates](#future-updates)

---

## Overview

The **Space Sciences and Astrodynamics** repository provides a robust toolkit for modeling, simulating, and analyzing space-based systems. It supports a wide range of functionalities including:
- Conversion between orbital elements and Cartesian state vectors.
- Simulation of spacecraft trajectories under gravitational and non-gravitational forces.
- Calculation of interplanetary trajectories and orbital perturbations.
- Numerical solutions to astrodynamics equations using advanced algorithms.
- Visualization of orbits, ground tracks, and other dynamic phenomena.

The toolkit is designed for educational and research purposes and is a resource for those interested in understanding or working within astrodynamics to explore and develop a deeper understanding of orbital mechanics.

## Requirements

### MATLAB
- **MATLAB**: Version R2021a or later recommended.
- Built-in functions such as `ode45` for numerical integration.
- `rsmooth.m` for proper functionality of dependent scripts.

### Python
- **Python**: Version 3.8 or later recommended
- Required packages:
  - `numpy`
  - `scipy`
  - `matplotlib`
- **Python**-translated **MATLAB** functions:
  - `atmosisa`
  - `quatinv`
  - `quatmultiply`
 
### C++
- **C++**: Version 23 or later recommended.
- **C++**-translated **MATLAB** functions:
  - `atmosisa`
  - `quatinv`
  - `quatmultiply`
  - `rsmooth`

## Features
- **Coordinate Transformations**: Convert between orbital elements and state vectors.
- **Orbit Propagation**: Simulate two-body and three-body systems.
- **Trajectory Analysis**: Compute interplanetary transfers and lunar trajectories.
- **Numerical Methods**: Implement advanced solvers like Runge-Kutta, Heun’s method, and bisection.
- **Ephemeris Calculations**: Predict positions of celestial bodies for a given Julian date.
- **Perturbation Modeling**: Analyze gravitational (e.g., J2 effects) and non-gravitational forces (e.g., solar radiation pressure).

## Files and Modules

### Course Code Modules

---

#### `Coordinate Transformation.py` 
- Converts between orbital elements and state vectors for orbital dynamics, calculating angular momentum, eccentricity, inclination, and other essential values.

---

#### `Family_of_orbits.py` 
- Simulates and visualizes a family of orbits with varying semi-major axis and eccentricity. Also includes an animation function for orbit visualization.

---

#### `Orbital Elements.py`
- Calculates orbital parameters from initial state vectors, such as vector magnitude, radial velocity, and angular momentum.

---

#### `Orbit_Simulator.py`
- Simulates orbital motion with transformations between orbital elements and Cartesian coordinates, providing a 3D orbit visualization.

<br>

### MATLAB Python and C++ Code Modules

---

#### `atan2d_0_360.m` | `atan2d_0_360.py` | `atan2d_0_360.cpp`

- The `atan2d_0_360` function is a specialized variant of the traditional `arctan` function, designed for scenarios where:
  - Angles need to be normalized to the range [0, 360] instead of the standard [-180, 180].
  - Consistent angle representation is required for astrodynamics, robotics, and other applications.
  - Improved readability and usability of angular calculations are beneficial in workflows.

- The test scripts included in both MATLAB and Python versions validate the `atan2d_0_360` function across a variety of input cases. These scripts:
  - Test the function with edge cases such as zero and axis-aligned vectors.
  - Evaluate performance on quadrants to ensure correct angle normalization.
  - Demonstrate the function's application in common scenarios.

![atan2d_0_360_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/atan2d_0_360_plot.png)

---

#### `atmosisa.py` | `atmosisa.cpp`

- The `atmosisa` function calculates atmospheric properties using an approximation of the International Standard Atmosphere (ISA) model. It provides:
  - Temperature (K) based on altitude in standard atmospheric layers.
  - Pressure (Pa) and density (kg/m³) for accurate environmental modeling.
  - Speed of sound (m/s) for applications in aerodynamics and astrodynamics.

- Key features include:
  - Support for altitudes from sea level up to approximately 85 km, with proper handling of stratospheric and mesosphere conditions.
  - Accurate representation of temperature gradients and pressure changes in each atmospheric layer.
  - User-friendly implementation that requires only altitude as input.

![atmosisa_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/atmosisa_plot.png)

---

#### `atmosphere.m` | `atmosphere.py` | `atmosphere.cpp`

- The `atmosphere` function calculates atmospheric density for altitudes ranging from sea level to 1000 km using exponential interpolation. This provides:
  - Accurate density values based on the U.S. Standard Atmosphere 1976 model.
  - Support for high-altitude atmospheric modeling for aerospace and astrodynamics applications.
  - Smooth transitions in density across altitude intervals for improved numerical stability.

- The test scripts included in both MATLAB and Python versions validate the `atmosphere` function and illustrate its functionality by:
  - Generating density profiles over the full altitude range.
  - Visualizing the relationship between altitude and density with logarithmic plots.
  - Demonstrating the smoothness and accuracy of the exponential interpolation.

![atmosphere_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/atmosphere_plot.png)

---

#### `bisect.m` | `bisect.py` | `bisect.cpp`

- The `bisect` function implements the bisection method for root-finding. This numerical algorithm is designed for scenarios where:
  - The root of a continuous function needs to be located within a specified interval.
  - High precision is required, with results computed to within a user-defined tolerance.
  - The method's simplicity and robustness are ideal for solving equations in astrodynamics and related fields.

- The test scripts included in both MATLAB and Python versions demonstrate the functionality of the `bisect` function by:
  - Calculating the three collinear Lagrange points ($L_1$, $L_2$, $L_3$) for the Earth-Moon system.
  - Validating the function's accuracy using predefined intervals and sample equations.
  - Outputting the results to show the convergence of the bisection method for practical problems.

![bisect_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/bisect_plot.png)

---

#### `coe_from_sv.m` | `coe_from_sv.py` | `coe_from_sv.cpp`

- The `coe_from_sv` function calculates classical orbital elements (COEs) from a given state vector using Algorithm 4.1. This is essential for:
  - Translating position and velocity vectors into orbital parameters.
  - Characterizing orbits in terms of eccentricity, inclination, and other key metrics.
  - Enabling comprehensive analysis and simulation in astrodynamics.

- The test scripts included in both MATLAB and Python versions demonstrate the use of `coe_from_sv` by:
  - Computing orbital elements for a sample state vector in a geocentric equatorial frame.
  - Validating the function's accuracy with known inputs and outputs.
  - Displaying additional orbital parameters such as the period for elliptical orbits.

![coe_from_sv_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/coe_from_sv_plot.png)

---

#### `cowell.m` | `cowell.py` | `cowell.cpp`

- The `cowell` function implements numerical integration using MATLAB's `ode45` or Python's `solve_ivp` to model satellite motion under atmospheric drag. It is designed for scenarios where:
  - Detailed trajectory predictions are needed, accounting for atmospheric drag effects.
  - Orbital evolution over extended periods must be analyzed.
  - Accurate modeling of satellite re-entry or decay due to drag is required.

- The scripts included in both MATLAB and Python versions integrates the equations of motion, incorporating:
  - Atmospheric density variations with altitude.
  - Drag force acting on the satellite.
  - Gravitational effects from Earth.

![cowell_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/cowell_plot.png)

---

#### `dcm_from_q.m` | `dcm_from_q.py` | `dcm_from_q.cpp`

- The `dcm_from_q` function calculates the direction cosine matrix (DCM) from a quaternion. This is essential for:
  - Translating quaternion-based orientation into matrix form for applications in aerospace and robotics.
  - Enabling efficient rotation transformations in simulations and real-time systems.
  - Maintaining numerical stability and accuracy in orientation calculations.

- The test scripts included in both MATLAB and Python versions validate the `dcm_from_q` function by:
  - Computing DCMs for known quaternions representing identity and 180-degree rotations about principal axes.
  - Ensuring the function produces accurate and consistent results across multiple input cases.
  - Demonstrating the utility of quaternion-to-DCM conversions in practical scenarios.

![dcm_from_q_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/dcm_from_q_plot.png)

---

#### `dcm_to_euler.m` | `dcm_to_euler.py` | `dcm_to_euler.cpp`

- The `dcm_to_euler` function extracts the classical Euler angles from a direction cosine matrix (DCM) using the sequence $R_3(\gamma) \cdot R_1(\beta) \cdot R_3(\alpha)$. This is useful for:
  - Converting rotational transformations into intuitive angle-based representations.
  - Analyzing and visualizing orientation in astrodynamics and robotics.
  - Supporting systems requiring Euler angle decomposition from matrix representations.

- The test scripts included in both MATLAB and Python versions validate the `dcm_to_euler` function by:
  - Testing various input matrices, including identity and specific rotations about principal axes.
  - Ensuring accurate conversion to Euler angles across diverse scenarios.
  - Demonstrating practical applications of DCM-to-Euler transformations.

![dcm_to_euler_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/dcm_to_euler_plot.png)

---

#### `dcm_to_ypr.m` | `dcm_to_ypr.py` | `dcm_to_ypr.cpp`

- The `dcm_to_ypr` function extracts the yaw, pitch, and roll angles from a direction cosine matrix (DCM) using the sequence $$R_1(\text{roll}) \cdot R_2(\text{pitch}) \cdot R_3(\text{yaw}) $$. This is useful for:
  - Converting matrix-based rotations into intuitive yaw, pitch, and roll angles.
  - Applications in aerospace, robotics, and navigation systems that require angle-based orientation.
  - Analyzing orientation data in a more comprehensible format.

- The test scripts included in both MATLAB and Python versions validate the `dcm_to_ypr` function by:
  - Testing various input DCMs, including identity and specific rotations about principal axes.
  - Ensuring accurate conversion to yaw, pitch, and roll angles across diverse cases.
  - Demonstrating the utility of DCM-to-YPR transformations in practical scenarios.

![dcm_to_ypr_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/dcm_to_ypr_plot.png)

---

#### `encke.m` | `encke.py` | `encke.cpp`

- The `encke` function uses Encke's method to numerically integrate the equations of motion for a satellite under the influence of $J_2$ gravitational perturbations. This code is designed for scenarios where:
  - Precise orbit propagation is required over long periods while accounting for perturbative forces.
  - Modeling the effects of Earth's oblateness $$( J_2 )$$ on orbital elements is essential.
  - Computational efficiency is important for handling small perturbations in satellite motion.

- The scripts included in both MATLAB and Python versions integrates the motion and computes variations in key orbital elements, such as:
  - Right ascension of the ascending node $$( \Delta\Omega )$$.
  - Argument of perigee $$( \Delta\omega )$$.
  - Eccentricity $$( \Delta e )$$, inclination $$( \Delta i )$$, and angular momentum $$( \Delta h )$$.

![encke_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/encke_plot1.png)

![encke_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/encke_plot2.png)

![encke_plot3.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/encke_plot3.png)

---

#### `f_and_g.m` | `f_and_g.py` | `f_and_g.cpp`

- The `f_and_g` function calculates the Lagrange coefficients $f$ and $g$, which are essential for solving orbital propagation problems. These coefficients:
  - Aid in determining the position and velocity of a celestial body after a time interval.
  - Incorporate the universal anomaly $x$ and other orbital parameters for precise calculations.
  - Support astrodynamics workflows requiring accurate orbital state predictions.

- The test scripts included in both MATLAB and Python versions validate the `f_and_g` function by:
  - Testing sample inputs for radial position, time, universal anomaly, and semimajor axis reciprocal.
  - Comparing calculated $f$ and $g$ coefficients to expected results.
  - Demonstrating the function's reliability in typical orbital mechanics scenarios.

![f_and_g_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/f_and_g_plot.png)

---

#### `f_and_g_ta.m` | `f_and_g_ta.py` | `f_and_g_ta.cpp`

- The `f_and_g_ta` function calculates the Lagrange $f$ and $g$ coefficients based on the change in true anomaly ($\Delta \theta$). This function is particularly useful for:
  - Propagating orbits using true anomaly changes without requiring a full numerical integration.
  - Efficiently determining position and velocity at a new orbital location.
  - Supporting analytical orbital mechanics workflows.

- The test scripts included in both MATLAB and Python versions validate the `f_and_g_ta` function by:
  - Using predefined initial position and velocity vectors, along with true anomaly changes, to compute $f$ and $g$.
  - Demonstrating the consistency of outputs for typical orbital mechanics scenarios.
  - Providing sample computations to verify accuracy and functionality.

![f_and_g_ta_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/f_and_g_ta_plot.png)

---

#### `fDot_and_gDot.m` | `fDot_and_gDot.py` | `fDot_and_gDot.cpp`

- The `fDot_and_gDot` function calculates the time derivatives of the Lagrange $f$ and $g$ coefficients, which are critical for:
  - Determining changes in position and velocity over time in orbital mechanics.
  - Supporting precise orbit propagation and trajectory predictions.
  - Analyzing dynamic changes in orbital elements using universal anomaly and radial positions.

- The test scripts included in both MATLAB and Python versions validate the `fDot_and_gDot` function by:
  - Providing sample inputs for radial positions, universal anomaly, and semimajor axis reciprocal.
  - Verifying the accuracy of computed $\dot{f}$ and $\dot{g}$ coefficients.
  - Demonstrating practical use cases for time-dependent orbital calculations.

![fDot_and_gDot_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/fDot_and_gDot_plot.png)

---

#### `fDot_and_gDot_ta.m` | `fDot_and_gDot_ta.py` | `fDot_and_gDot_ta.cpp`

- The `fDot_and_gDot_ta` function calculates the time derivatives of the Lagrange $f$ and $g$ coefficients based on a change in the true anomaly ($\Delta \theta$). This is particularly useful for:
  - Analyzing time-dependent changes in orbital position and velocity.
  - Supporting orbital mechanics calculations involving angular momentum and radial velocity.
  - Providing efficient computations for trajectory propagation in astrodynamics.

- The test scripts included in both MATLAB and Python versions validate the `fDot_and_gDot_ta` function by:
  - Using sample initial position and velocity vectors and changes in true anomaly to compute $\dot{f}$ and $\dot{g}$.
  - Verifying the accuracy of outputs for practical orbital mechanics scenarios.
  - Demonstrating the function's reliability for analyzing orbital element variations.

![fDot_and_gDot_ta_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/fDot_and_gDot_ta_plot.png)

---

#### `gauss.m` | `gauss.py` | `gauss.cpp`

- The `gauss` function implements Gauss's method for preliminary orbit determination, including iterative improvement. This algorithm calculates the state vector (position and velocity) of an orbiting body from angles-only observations at three closely spaced times. Key features include:
  - Precise orbit determination using minimal observation data.
  - Iterative refinement to improve the accuracy of the calculated state vector.
  - Support for real-world applications in satellite tracking and orbital mechanics.

- The test scripts included in both MATLAB and Python versions validate the `gauss` function by:
  - Using sample observational data for right ascension, declination, and sidereal time.
  - Computing orbital elements such as angular momentum, eccentricity, inclination, and semimajor axis.
  - Comparing results with and without iterative improvement to demonstrate the method's convergence and reliability.

![gauss_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/gauss_plot.png)

---

#### `gibbs.m` | `gibbs.py` | `gibbs.cpp`

- The `gibbs` function implements Gibbs's method for preliminary orbit determination using three coplanar geocentric position vectors. It is designed for:
  - Calculating the velocity corresponding to the second position vector.
  - Determining orbital elements such as angular momentum, eccentricity, and inclination.
  - Supporting applications in satellite tracking and celestial mechanics.

- The test scripts included in both MATLAB and Python versions validate the `gibbs` function by:
  - Using sample coplanar position vectors to compute the velocity vector and orbital elements.
  - Checking for coplanarity of input vectors to ensure method applicability.
  - Demonstrating the accuracy and reliability of the function in typical orbital scenarios.

![gibbs_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/gibbs_plot.png)

---

#### `gravity_turn.m` | `gravity_turn.py` | `gravity_turn.cpp`

- The `gravity_turn` function numerically integrates the equations of motion for a gravity turn trajectory. This is particularly useful for:
  - Simulating launch vehicle trajectories under the influence of thrust, drag, and gravity.
  - Calculating dynamic parameters such as altitude, velocity, and flight path angle.
  - Analyzing velocity losses due to drag and gravity.

- The scripts included in both MATLAB and Python versions integrates the gravity turn trajectory, incorporating:
  - Dynamic simulation of pitchover and powered ascent phases.
  - Visualization of trajectory and dynamic pressure profiles.
  - Support for varying vehicle and atmospheric parameters to model real-world scenarios.

![gravity_turn_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/gravity_turn_plot1.png)

![gravity_turn_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/gravity_turn_plot2.png)

---

#### `ground_track.m` | `ground_track.py` | `ground_track.cpp`

- The `ground_track` function calculates and plots the ground track of a satellite using its orbital elements. This is useful for:
  - Visualizing the satellite's trajectory relative to Earth's surface over multiple orbits.
  - Simulating orbital evolution, including nodal regression and perigee precession due to $J_2$ perturbations.
  - Analyzing the satellite's coverage and ground visibility.

- The scripts included in both MATLAB and Python versions integrates the ground track, incorporating:
  - Computation of right ascension and declination histories for a specified time interval.
  - Separation of the ground track into distinct curves to account for orbital discontinuities.
  - Detailed plots of the ground track on Earth's surface for visualization.

![ground_track_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/ground_track_plot.png)

---

#### `heun.m` | `heun.py` | `heun.cpp`

- The `heun` function implements Heun's method, a predictor-corrector numerical integration technique, for solving systems of first-order ordinary differential equations (ODEs). This function is designed for:
  - Accurate numerical integration of time-dependent ODEs in engineering and physics applications.
  - Handling systems with stringent error tolerance and iterative correction.
  - Supporting user-defined differential equations through callable functions.

- The test scripts included in both MATLAB and Python versions demonstrate the use of `heun` by:
  - Solving a second-order damped single-degree-of-freedom spring-mass system under sinusoidal forcing.
  - Comparing solutions obtained with different time step sizes to illustrate convergence and accuracy.
  - Visualizing displacement responses over time for varying numerical resolutions.

![heun_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/heun_plot.png)

---

#### `integrate_thrust.m` | `integrate_thrust.py` | `integrate_thrust.cpp`

- The `integrate_thrust` function uses numerical integration to model the dynamics of a spacecraft during a delta-v burn. It calculates the resulting trajectory and orbital parameters after the burn. This function is ideal for:
  - Simulating orbital maneuvers, including changes in velocity and trajectory.
  - Computing the apogee and other orbital elements of the post-burn orbit.
  - Evaluating thrust, mass loss, and other dynamic factors affecting spacecraft motion.

- The scripts included in both MATLAB and Python versions integrates the thrust of the spacecraft, incorporating:
  - Accurate integration of motion equations using the Runge-Kutta-Fehlberg (RKF45) method.
  - Determination of state vectors and orbital elements before and after the burn.
  - Support for flexible inputs to simulate a variety of propulsion scenarios.

![integrate_thrust_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/integrate_thrust_plot.png)

---

#### `interplanetary.m` | `interplanetary.py` | `interplanetary.cpp`

- The `interplanetary` function calculates the trajectory of a spacecraft traveling from one planet to another using patched conic approximations. This function is ideal for:
  - Planning interplanetary missions by determining heliocentric state vectors at departure and arrival.
  - Computing the spacecraft's hyperbolic excess velocities at both planets.
  - Supporting mission design and analysis with key orbital elements and time-of-flight calculations.

- The scripts included in both MATLAB and Python versions demonstrate the use of `interplanetary` by:
  - Solving for trajectories between specified departure and arrival dates for Earth and Mars, and other planents in the Solar System.
  - Validating computed orbital elements and hyperbolic excess velocities.
  - Providing a comprehensive analysis of the interplanetary transfer.

![interplanetary_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/interplanetary_plot.png)

---

#### `J0.m` | `J0.py` | `J0.cpp`

- The `J0` function calculates the Julian day number at 0 UT (Universal Time) for a given date. This is essential for:
  - Converting calendar dates into a numerical format suitable for astronomical and orbital calculations.
  - Supporting time-dependent computations in astrodynamics and celestial mechanics.
  - Ensuring compatibility with standard equations for Julian date arithmetic.

- The test scripts included in both MATLAB and Python versions validate the `J0` function by:
  - Computing the Julian day number for a specific date and time.
  - Verifying the output against known examples to ensure accuracy.
  - Demonstrating the use of Julian day numbers in time-sensitive computations.

![J0_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/J0_plot.png)

---

#### `J2_perturbation.m` | `J2_perturbation.py` | `J2_perturbation.cpp`

- The `J2_perturbation` function numerically integrates Gauss's planetary equations to analyze the effects of Earth's oblateness ($J_2 $) on the orbital elements of a satellite. This function is designed for:
  - Investigating long-term perturbations in orbital parameters caused by $J_2$.
  - Simulating changes in right ascension, inclination, eccentricity, and argument of perigee.
  - Supporting mission planning and satellite lifetime analysis with precise perturbation modeling.

- The scripts included in both MATLAB and Python versions include:
  - Integration of the variational equations for detailed time histories of orbital elements.
  - Visualization of the effects of $J_2$ on angular momentum, eccentricity, and inclination.
  - Flexible setup for analyzing different orbital regimes and initial conditions.

![J2_perturbation_plot_1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/J2_perturbation_plot_1.png)

![J2_perturbation_plot_2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/J2_perturbation_plot_2.png)

![J2_perturbation_plot_3.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/J2_perturbation_plot_3.png)

---

#### `kepler_E.m` | `kepler_E.py` | `kepler_E.cpp`

- The `kepler_E` function solves Kepler's equation $E - e \sin(E) = M$ for the eccentric anomaly ($E$) using Newton's method. This function is essential for:
  - Determining orbital positions in elliptical trajectories.
  - Supporting astrodynamics calculations where precise orbital element computation is required.
  - Converting mean anomaly ($M$) to eccentric anomaly ($E$) for orbital mechanics analyses.

- The test scripts included in both MATLAB and Python versions validate the `kepler_E` function by:
  - Solving Kepler's equation for a given eccentricity and mean anomaly.
  - Verifying the convergence of Newton's method to within a specified tolerance.
  - Demonstrating the practical application of the function in orbital simulations.

![kepler_E_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/kepler_E_plot.png)

---

#### `kepler_H.m` | `kepler_H.py` | `kepler_H.cpp`

- The `kepler_H` function solves Kepler's equation for hyperbolic orbits, $e \sinh(F) - F = M$, using Newton's method to compute the hyperbolic eccentric anomaly ($F$). This function is essential for:
  - Determining orbital positions in hyperbolic trajectories.
  - Supporting interplanetary mission design involving escape or flyby trajectories.
  - Converting hyperbolic mean anomaly ($M$) to hyperbolic eccentric anomaly ($F$).

- The test scripts included in both MATLAB and Python versions validate the `kepler_H` function by:
  - Solving Kepler's equation for a given eccentricity and hyperbolic mean anomaly.
  - Ensuring the convergence of Newton's method to the desired accuracy.
  - Demonstrating the utility of the function in analyzing hyperbolic orbits.

![kepler_H_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/kepler_H_plot.png)

---

#### `kepler_U.m` | `kepler_U.py` | `kepler_U.cpp`

- The `kepler_U` function solves the universal Kepler's equation to compute the universal anomaly ($x$) using Newton's method. This function is crucial for:
  - Propagating orbits in any conic section (elliptical, parabolic, or hyperbolic).
  - Calculating position and velocity vectors over time for celestial bodies or spacecraft.
  - Enabling generalized solutions for orbital mechanics problems independent of orbit type.

- The test scripts included in both MATLAB and Python versions validate the `kepler_U` function by:
  - Solving the universal Kepler's equation for a given radial position, velocity, and elapsed time.
  - Verifying convergence and accuracy of the computed universal anomaly.
  - Demonstrating the function's application in dynamic orbit propagation scenarios.

![kepler_U_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/kepler_U_plot.png)

---

#### `lambert.m` | `lambert.py` | `lambert.cpp`

- The `lambert` function solves Lambert's problem, which calculates the velocity vectors required to transfer a spacecraft between two position vectors in a specified time. This is essential for:
  - Planning interplanetary and orbital transfer trajectories.
  - Computing initial and final velocities for given orbital configurations.
  - Supporting mission design and analysis in astrodynamics.

- The test scripts included in both MATLAB and Python versions validate the `lambert` function by:
  - Solving for velocity vectors for given initial and final position vectors and time of flight.
  - Calculating and verifying orbital elements, including angular momentum, eccentricity, and inclination.
  - Demonstrating the application of the function in single-revolution transfer orbits.

![lambert_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/lambert_plot.png)

---

#### `los.m` | `los.py` | `los.cpp`

- The `los` function determines whether a satellite is in Earth's shadow based on the satellite's and Sun's Earth-Centered Inertial (ECI) position vectors. This function is critical for:
  - Assessing solar exposure for satellites in orbit.
  - Determining periods of Earth shadow for satellite operations.
  - Supporting power and thermal analysis for spacecraft.

- The test scripts included in both MATLAB and Python versions validate the `los` function by:
  - Simulating scenarios where the satellite is in sunlight and Earth's shadow.
  - Analyzing different satellite positions, including edge cases such as satellites near Earth's surface.
  - Verifying function accuracy in determining the line-of-sight status between the satellite and Sun.

![los_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/los_plot.png)

---

#### `LST.m` | `LST.py` | `LST.cpp`

- The `LST` function calculates the Local Sidereal Time (LST) at a specified location on Earth. This is essential for:
  - Determining the orientation of the Earth relative to celestial objects.
  - Supporting astronomical observations and satellite tracking.
  - Converting time and geographic location into a sidereal time format for astrodynamics computations.

- The test scripts included in both MATLAB and Python versions validate the `LST` function by:
  - Computing LST for a given date, Universal Time (UT), and geographic location.
  - Echoing input data and comparing computed LST values to expected results.
  - Demonstrating the use of LST for observing celestial events or satellite trajectories.

![LST_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/LST_plot.png)

---

#### `lunar_perturbation_test.m` | `lunar_perturbation_test.py` | `lunar_perturbation_test.cpp`

- The `lunar_perturbation_test` function models the effects of the Moon's gravitational perturbations on satellite orbits using the Gauss variational equations. This function is ideal for:
  - Analyzing long-term orbital changes due to lunar gravity.
  - Simulating orbital element variations such as inclination, right ascension, and argument of perigee.
  - Supporting mission planning and stability assessments for Earth-orbiting satellites.

- The test scripts included in both MATLAB and Python versions validate the lunar perturbations, by:
  - Numerical integration of Gauss's variational equations using MATLAB's `ode45` and Python's `solve_ivp`.
  - Analysis of three orbital regimes: geostationary (GEO), highly elliptical (HEO), and low Earth orbit (LEO).
  - Visualization of orbital element evolution under lunar perturbations.

![lunar_perturbation_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/lunar_perturbation_plot1.png)

![lunar_perturbation_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/lunar_perturbation_plot2.png)

---

#### `lunar_position.m` | `lunar_position.py` | `lunar_position.cpp`

- The `lunar_position` function calculates the geocentric equatorial position vector of the Moon for a given Julian date. This is essential for:
  - Determining the Moon's position relative to Earth in celestial mechanics.
  - Supporting lunar mission planning and astrodynamics calculations.
  - Providing precise Moon position data for time-sensitive orbital computations.

- The test scripts included in both MATLAB and Python versions validate the `lunar_position` function by:
  - Calculating the Moon's position vector for key Julian dates, including J2000 and specific historical dates.
  - Comparing the results with known values to ensure accuracy.
  - Demonstrating the function's utility in lunar and Earth-centered orbital analyses.

![lunar_position_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/lunar_position_plot.png)

---

#### `lunar_trajectory.m` | `lunar_trajectory.py` | `lunar_trajectory.cpp`

- The `lunar_trajectory` function models the trajectory of a spacecraft under the gravitational influence of Earth and the Moon using a restricted three-body problem framework. Key features include:
  - Accurate modeling of spacecraft motion with initial conditions and parameters.
  - Visualization of the trajectory in both the Earth-Centered Inertial (ECI) frame and a Moon-fixed frame.
  - Calculation of key points, such as perilune and lunar arrival conditions, for trajectory analysis.
  - Support for mission planning, trajectory optimization, and astrodynamics research.

- The test scripts included in both MATLAB and Python versions leverage `ode45` in MATLAB and `solve_ivp` in Python for numerical integration. The functions include detailed graphical outputs and calculated parameters, such as:
  - Moon's position and velocity.
  - Probe's perilune altitude and relative speed.
  - Trajectory inclination and osculating plane details.

![lunar_trajectory_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/lunar_trajectory_plot1.png)

![lunar_trajectory_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/lunar_trajectory_plot2.png)

![lunar_trajectory_plot3.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/lunar_trajectory_plot3.png)

---

#### `month_planet_names.m` | `month_planet_names.py` | `month_planet_names.cpp`

- The `month_planet_names` function returns the name of a month and a planet based on their respective numeric IDs. This utility is designed for:
  - Mapping numeric IDs to human-readable month and planet names.
  - Supporting educational tools or software requiring quick access to such mappings.
  - Simplifying conversions between numerical representations and their corresponding names.

- The test scripts included in both MATLAB and Python versions validate the `month_planet_names` function by:
  - Iterating through all valid month and planet IDs.
  - Printing the corresponding names for verification.
  - Demonstrating accurate mappings and robustness across valid inputs.

![month_planet_names_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/month_planet_names_plot.png)

---

#### `orbit.m` | `orbit.py` | `orbit.cpp`

- The `orbit` function computes the orbit of a spacecraft by numerically solving the two-body problem using the Runge-Kutta-Fehlberg (RKF45) method. This function is designed for:
  - Simulating orbital trajectories and calculating spacecraft positions and velocities over time.
  - Determining key orbital parameters, including maximum and minimum radii and the corresponding velocities.

- The scripts included in both MATLAB and Python versions validate the `orbit` function by:
  - Input flexibility for initial conditions, such as position and velocity vectors.
  - Numerical integration of the two-body equations of motion.
  - Detailed plots of the orbital trajectory, including annotations for critical points such as the maximum and minimum altitudes.

![orbit_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/orbit_plot.png)

---

#### `planet_elements_and_sv.m` | `planet_elements_and_sv.py` | `planet_elements_and_sv.cpp`

- The `planet_elements_and_sv` function calculates the heliocentric orbital elements and state vectors (position and velocity) of a planet for a given date and time. This function is essential for:
  - Determining the current position and motion of planets in the solar system.
  - Supporting astrodynamics and celestial mechanics computations.
  - Providing data for mission planning and planetary observation.

- The test scripts included in both MATLAB and Python versions validate the `planet_elements_and_sv` function by:
  - Calculating orbital elements and state vectors for planets such as Earth and Mars.
  - Verifying Julian date calculations and their impact on planetary positions.
  - Demonstrating the accuracy of heliocentric position and velocity computations.

![planet_elements_and_sv_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/planet_elements_and_sv_plot.png)

---

#### `q_from_dcm.m` | `q_from_dcm.py` | `q_from_dcm.cpp`

- The `q_from_dcm` function calculates the quaternion representation of a given direction cosine matrix (DCM). This is essential for:
  - Converting rotation matrices into quaternions for use in spacecraft attitude dynamics and robotics.
  - Enabling efficient computation of orientation changes.
  - Supporting applications in control systems, simulation, and 3D transformations.

- The test scripts included in both MATLAB and Python versions validate the `q_from_dcm` function by:
  - Testing with known DCMs, such as the identity matrix and 180-degree rotations about the principal axes.
  - Printing and comparing the computed quaternions to expected results.
  - Demonstrating the function's accuracy and reliability for practical rotation matrix conversions.

![q_from_dcm_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/q_from_dcm_plot.png)

---

#### `quat_rotate.m` | `quat_rotate.py` | `quat_rotate.cpp`

- The `quat_rotate` function rotates a vector by a unit quaternion, providing a method for:
  - Transforming vectors in 3D space using quaternion-based rotation.
  - Supporting spacecraft attitude simulations and 3D graphics applications.
  - Simplifying rotational transformations with efficient quaternion mathematics.

- The test scripts included in both MATLAB and Python versions validate the `quat_rotate` function by:
  - Rotating vectors with predefined quaternions representing rotations around principal axes.
  - Comparing rotated vectors against expected results to ensure accuracy.
  - Demonstrating the function's application in 3D rotational computations.

#### `quatinv.py` and `quatmultiply.py` | `quatinv.cpp` and `quatmultiply.cpp`

- These Python scripts implement the functionality of MATLAB's built-in `quatinv` and `quatmultiply` functions. They are essential for:
  - Calculating the inverse of a quaternion with `quatinv`.
  - Performing quaternion multiplication with `quatmultiply`.

- These scripts replicate MATLAB's behavior, enabling seamless transition for MATLAB users moving to Python environments. They are utilized in functions such as `quat_rotate` for quaternion-based calculations.

![quat_rotate_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/quat_rotate_plot.png)

---

#### `ra_and_dec_from_r.m` | `ra_and_dec_from_r.py` | `ra_and_dec_from_r.cpp`

- The `ra_and_dec_from_r` function calculates the right ascension and declination from a geocentric equatorial position vector. This function is essential for:
  - Converting position vectors into angular coordinates for use in celestial mechanics and navigation.
  - Determining the orientation of an object in Earth's sky from its position in space.
  - Supporting astrodynamics applications, including satellite tracking and orbital visualization.

- The test scripts included in both MATLAB and Python versions validate the `ra_and_dec_from_r` function by:
  - Computing right ascension and declination for a sample position vector.
  - Comparing the results to expected values to ensure accuracy.
  - Demonstrating the practical application of the function in orbital mechanics.

![ra_and_dec_from_r_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/ra_and_dec_from_r_plot.png)

---

#### `relative_motion.m` | `relative_motion.py` | `relative_motion.cpp`

- The `relative_motion` function calculates and plots the motion of a chaser satellite (B) relative to a target satellite (A) in orbit. This function is useful for:
  - Simulating relative motion dynamics between two satellites for formation flying or rendezvous missions.
  - Visualizing relative trajectories in a co-moving reference frame.
  - Supporting analysis of proximity operations in astrodynamics.

- The scripts included in both MATLAB and Python versions include:
  - Numerical integration of relative motion equations using MATLAB's `ode45` or Python's `rkf45`.
  - Representation of the relative trajectory of the chaser in the co-moving frame of the target.
  - Incorporation of orbital parameters, including angular momentum, eccentricity, and mean motion, to accurately model the relative dynamics.

![relative_motion_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/relative_motion_plot.png)

---

#### `rk1_4.m` | `rk1_4.py` | `rk1_4.cpp`

- The `rk1_4` function implements the Runge-Kutta numerical integration methods (RK1, RK2, RK3, and RK4) to solve systems of first-order ordinary differential equations $\frac{dy}{dt} = f(t, y)$. This function is essential for:
  - Performing numerical simulations of dynamic systems.
  - Providing flexibility to choose between integration methods depending on accuracy and computational efficiency.
  - Supporting the analysis of time-dependent problems in astrodynamics and engineering.

- The test scripts included in both MATLAB and Python versions validate the `rk1_4` function by:
  - Solving a damped spring-mass system with sinusoidal forcing for different time steps and integration orders.
  - Comparing numerical solutions to the exact analytical response.
  - Demonstrating the relative accuracy and stability of the four Runge-Kutta methods.

![rk1_4_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/rk1_4_plot.png)

---

#### `rkf45.m` | `rkf45.py` | `rkf45.cpp`

- The `rkf45` function implements the Runge-Kutta-Fehlberg 4(5) method, a numerical integration algorithm with adaptive step size control. It is designed for:
  - Solving systems of first-order ordinary differential equations $\frac{dy}{dt} = f(t, y)$ with high accuracy and efficiency.
  - Applications in astrodynamics, engineering, and physics requiring precise integration over varying time steps.
  - Supporting problems with stringent error tolerance and dynamic step size adjustments.

- The test scripts included in both MATLAB and Python versions validate the `rkf45` function by:
  - Solving the two-body rectilinear motion and restricted three-body problem.
  - Demonstrating the function's application to satellite motion and spacecraft trajectory analysis.
  - Comparing numerical results to analytical expectations and ensuring adaptive step control.

![rkf45_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/rkf45_plot.png)

![rkf45_lunar_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/rkf45_lunar_plot.png)

---

#### `rv_from_observe.m` | `rv_from_observe.py` | `rv_from_observe.cpp`

- The `rv_from_observe` function calculates the geocentric equatorial position and velocity vectors of an object from radar observations, including range, azimuth, elevation angle, and their rates. This function is essential for:
  - Translating radar observation data into state vectors for orbit determination.
  - Supporting satellite tracking and orbit prediction applications.
  - Providing a foundation for converting observational data into orbital elements.

- The test scripts included in both MATLAB and Python versions validate the `rv_from_observe` function by:
  - Using simulated radar observations to compute state vectors.
  - Applying the state vectors to derive orbital elements such as angular momentum, eccentricity, and inclination.
  - Comparing results to known values to ensure correctness and demonstrating the function's utility in orbit determination tasks.

![rv_from_observe_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/rv_from_observe_plot.png)

---

#### `rv_from_r0v0.m` | `rv_from_r0v0.py` | `rv_from_r0v0.cpp`

- The `rv_from_r0v0` function computes the state vector (position and velocity) of a spacecraft at a future time, given its initial state vector and the elapsed time. This is fundamental for:
  - Propagating orbital trajectories using Kepler's equations and universal anomaly methods.
  - Simulating spacecraft motion in the two-body problem framework.
  - Supporting orbital dynamics analysis and mission planning tasks.

- The test scripts included in both MATLAB and Python versions validate the `rv_from_r0v0` function by:
  - Propagating a predefined initial state vector through a given time interval.
  - Calculating and echoing the final state vector for accuracy validation.
  - Demonstrating the function's utility in dynamic orbit propagation tasks.

![rv_from_r0v0_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/rv_from_r0v0_plot.png)

---

#### `rv_from_r0v0_ta.m` | `rv_from_r0v0_ta.py` | `rv_from_r0v0_ta.cpp`

- The `rv_from_r0v0_ta` function calculates the state vector (position and velocity) of a spacecraft after a specified change in true anomaly, given the initial state vector. This function is essential for:
  - Propagating orbital trajectories based on changes in true anomaly.
  - Determining spacecraft position and velocity at specific points in the orbit.
  - Supporting astrodynamics analysis and mission design tasks.

- The test scripts included in both MATLAB and Python versions validate the `rv_from_r0v0_ta` function by:
  - Using predefined initial conditions and changes in true anomaly to compute final state vectors.
  - Comparing the computed results to expected values for accuracy verification.
  - Demonstrating the function's utility in practical orbital mechanics scenarios.

![rv_from_r0v0_ta_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/rv_from_r0v0_ta_plot.png)

---

#### `rva_relative.m` | `rva_relative.py` | `rva_relative.cpp`

- The `rva_relative` function calculates the position, velocity, and acceleration of a spacecraft B relative to spacecraft A in A's Local Vertical Local Horizontal (LVLH) frame. This function is vital for:
  - Analyzing relative motion in rendezvous and proximity operations.
  - Simulating formation flying and cooperative satellite missions.
  - Supporting mission planning for docking and station-keeping tasks.

- The test scripts included in both MATLAB and Python versions validate the `rva_relative` function by:
  - Using orbital parameters to compute the relative position, velocity, and acceleration of one spacecraft with respect to another.
  - Comparing outputs to expected values derived from analytical solutions.
  - Demonstrating its application in dynamic simulation scenarios.

![rva_relative_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/rva_relative_plot.png)

---

#### `simpsons_lunar_ephemeris.m` | `simpsons_lunar_ephemeris.py` | `simpsons_lunar_ephemeris.cpp`

- The `simpsons_lunar_ephemeris` function computes the state vector (position and velocity) of the Moon relative to Earth's geocentric equatorial frame for a given Julian date. This function is based on a curve fit to JPL's DE200 ephemeris model, providing:
  - Accurate lunar position and velocity data for astrodynamics applications.
  - Essential inputs for spacecraft mission planning and lunar navigation.
  - A simplified, efficient approach suitable for onboard flight software.

- The test scripts included in both MATLAB and Python versions validate the `simpsons_lunar_ephemeris` function by:
  - Computing lunar state vectors for a sample Julian date and verifying the results.
  - Displaying the position and velocity vectors for inspection.
  - Demonstrating the function's application in lunar ephemeris modeling.

![simpsons_lunar_ephemeris_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/simpsons_lunar_ephemeris_plot1.png)

---

#### `solar_perturbation.m` | `solar_perturbation.py` | `solar_perturbation.cpp`

- The `solar_perturbation` function simulates the effects of solar gravitational perturbation on a spacecraft's orbit by integrating the Gauss variational equations. This function is designed for:
  - Analyzing the long-term influence of solar gravity on satellite orbits.
  - Supporting mission planning by quantifying perturbations on different orbital regimes (e.g., LEO, HEO, GEO).
  - Providing insights into orbital element variations due to external gravitational effects.

- The scripts included in both MATLAB and Python versions generate comparative simulations for multiple orbits. The function:
  - Computes changes in right ascension, inclination, and argument of perigee over time.
  - Utilizes advanced integration techniques to achieve high accuracy.
  - Visualizes the impact of solar perturbations on orbital elements through detailed plots.

![solar_perturbation_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/solar_perturbation_plot1.png)

![solar_perturbation_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/solar_perturbation_plot2.png)

---

#### `solar_position.m` | `solar_position.py` | `solar_position.cpp`

- The `solar_position` function calculates the geocentric equatorial position vector of the Sun for a given Julian date. This is crucial for:
  - Determining the Sun's position relative to Earth in celestial mechanics and satellite mission planning.
  - Supporting solar radiation pressure calculations and perturbation modeling.
  - Providing accurate inputs for astrodynamics and space mission design.

- The test scripts included in both MATLAB and Python versions validate the `solar_position` function by:
  - Calculating the Sun's position for specific Julian dates, such as the Winter Solstice.
  - Displaying results, including the apparent ecliptic longitude, obliquity of the ecliptic, and geocentric position vector.
  - Demonstrating the function's utility for solar ephemeris computations.

![solar_position_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/solar_position_plot.png)

---

#### `solar_radiation_pressure.m` | `solar_radiation_pressure.py` | `solar_radiation_pressure.cpp`

- The `solar_radiation_pressure` function models the effects of solar radiation pressure on a satellite's orbital elements by solving the Gauss planetary equations. This is particularly useful for:
  - Evaluating the long-term influence of solar radiation on satellite orbits.
  - Supporting mission planning by quantifying non-gravitational perturbations.
  - Analyzing changes in angular momentum, eccentricity, and other orbital parameters over time.

- The scripts included in both MATLAB and Python versions provide capabilities to:
  - Simulate variations in orbital elements caused by solar pressure.
  - Visualize results through plots of key orbital parameters.
  - Highlight differences between initial and perturbed orbital states.

![solar_radiation_pressure_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/solar_radiation_pressure_plot.png)

![solar_radiation_pressure_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/solar_radiation_pressure_plot2.png)

---

#### `spinning_top.m` | `spinning_top.py` | `spinning_top.cpp`

- The `spinning_top` function numerically integrates Euler's equations of motion for a spinning top. It uses quaternions to track the time evolution of the top's orientation. This implementation is designed for:
  - Understanding the rotational dynamics of a rigid body under gravity.
  - Analyzing Euler angles, angular velocities, and the effects of moments of inertia.
  - Providing insights into the nutation, precession, and spin behavior of spinning tops.

- The scripts included in both MATLAB and Python versions provide:
  - A quaternion-based approach for orientation tracking.
  - RKF45 numerical ODE solver for integrating the equations of motion.
  - Auxiliary functions for converting between direction cosine matrices (DCM) and Euler angles or quaternions.
  - Precession, nutation, and spin angles over time.
  - Corresponding angular rates and visualizations to highlight dynamic changes.

![spinning_top_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/spinning_top_plot1.png)

![spinning_top_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/spinning_top_plot2.png)

![spinning_top_plot3.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/spinning_top_plot3.png)

---

#### `stumpC.m` | `stumpC.py` | `stumpC.cpp`

- The `stumpC` function evaluates the Stumpff function $C(z)$, a key mathematical function in orbital mechanics, particularly for solving Kepler's equation using the universal variable formulation. This implementation is designed for:
  - Efficiently handling cases for positive, negative, and zero values of $z$.
  - Supporting numerical algorithms in astrodynamics that require Stumpff function computations.
  - Providing robust and reliable results for diverse orbital scenarios.

- The test scripts included in both MATLAB and Python versions validate the `stumpC` function by:
  - Testing it with a range of positive, negative, and zero values of $z$.
  - Verifying known outcomes for specific cases to ensure accuracy.
  - Displaying results in a structured format for interpretation and comparison.

![stumpC_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/stumpC_plot.png)

---

#### `stumpS.m` | `stumpS.py` | `stumpS.cpp`

- The `stumpS` function computes the Stumpff function $S(z)$, essential in solving Kepler's equation using the universal variable formulation. It handles different cases of $z$ (positive, negative, and zero) to provide:
  - Accurate evaluations for orbital mechanics problems.
  - Support for numerical algorithms in astrodynamics requiring Stumpff function computations.
  - Reliable results across various orbital scenarios.

- The test scripts included in both MATLAB and Python versions validate the `stumpS` function by:
  - Testing it with a range of positive, negative, and zero values of $z$.
  - Comparing results against known outcomes for specific values.
  - Ensuring robust and accurate implementation for all cases.

![stumpS_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/stumpS_plot.png)

---

#### `sv_from_coe.m` | `sv_from_coe.py` | `sv_from_coe.cpp`

- The `sv_from_coe` function computes the state vector (position and velocity) from classical orbital elements. This is particularly useful in astrodynamics for:
  - Translating orbital elements like eccentricity, inclination, and true anomaly into precise position and velocity vectors in a geocentric equatorial frame.
  - Supporting simulations and analyses of orbital mechanics, mission planning, and trajectory design.
  - Ensuring consistent and accurate transformations between orbital and Cartesian coordinates.

- The test scripts included in both MATLAB and Python versions validate the `sv_from_coe` function across a variety of cases:
  - Test the function with predefined orbital element sets to ensure accuracy.
  - Evaluate the generated state vectors for correctness against theoretical expectations.

![sv_from_coe_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/sv_from_coe_plot.png)

---

#### `threebody.m` | `threebody.py` | `threebody.cpp`

- The `threebody` function provides a graphical simulation of the motion of three celestial bodies interacting under gravitational forces in a two-dimensional plane. It uses:
  - Newton's laws of motion and the gravitational constant to compute positions and velocities over time.
  - Numerical integration (Runge-Kutta in MATLAB and Euler's method in Python) to solve the system of differential equations governing the bodies' motion.

- The scripts included in both MATLAB and Python versions validates:
  - The center of mass and relative positions of the bodies.
  - The trajectories of the bodies in an inertial frame and relative to the center of mass.
  - Useful for studying three-body problems in celestial mechanics, such as planetary motion or gravitational dynamics in a triple star system.

![threebody_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/threebody_plot1.png)

![threebody_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/threebody_plot2.png)

---

#### `twobody3d.py` | `twobody3d.py` | `twobody3d.cpp`

- The `twobody3d` function numerically solves the two-body problem in three dimensions relative to an inertial frame. The simulation uses the RKF4(5) numerical integration method to model the gravitational interactions between two massive bodies.
  - Inputs include initial positions and velocities of the bodies, their masses, and the time interval of the simulation.
  - Outputs provide trajectories, velocities, and the center of mass for the system at discrete time intervals.

- The scripts included in both MATLAB and Python versions validates:
  - The motion of each body relative to an inertial frame.
  - The motion of the second body and the center of mass relative to the first body.
  - The motion of both bodies relative to the center of mass.

![twobody3d_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/twobody3d_plot1.png)

![twobody3d_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/twobody3d_plot2.png)

![twobody3d_plo3t.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/Plots/twobody3d_plot3.png)

---

## Usage

### MATLAB
Clone the repository and navigate to the MATLAB files directory. Navigate to the .m files and execute scripts as needed.

### Python
Clone the repository and navigate to the Python files directory. Install the required .py packages and run scripts directly.

### C++
Clone the repository and navigate to the C++ and C++ Header files directory. Navigate to the .cpp and .h files and execute scripts as needed.

## Installation

Clone the repository:
```bash
git clone https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics.git
```

---

## Future Updates
Actively working on extending functionality to other programming languages, including:
- **Fortran**: Integration with legacy systems and specialized scientific computations.
- **Plot Scripts**: Will add the Python Plot Scripts in the near future for replicating/modifying plots displayed in README.

Additionally, ongoing improvements will be made to enhance existing modules, improve usability, and incorporate new features as needed.
