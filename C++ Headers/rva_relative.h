// ALGORITHM 7.1: FIND THE POSITION, VELOCITY, AND ACCELERATION
// OF B RELATIVE TO A'S LVLH FRAME
#ifndef RVA_RELATIVE_H
#define RVA_RELATIVE_H

#include <iostream>
#include <vector>
#include <cmath>
#include "sv_from_coe.h"

/*
    This function uses the state vectors of spacecraft A and B
    to find the position, velocity and acceleration of B relative
    to A in the LVLH frame attached to A (see Figure 7.1).

    rA,vA       - state vector of A (km, km/s)
    rB,vB       - state vector of B (km, km/s)
    mu          - gravitational parameter (km^3/s^2)
    hA          - angular momentum vector of A (km^2/s)
    i, j, k     - unit vectors along the x, y and z axes of A's
                  LVLH frame
    QXx         - DCM of the LVLH frame relative to the geocentric
                  equatorial frame (GEF)
    Omega       - angular velocity of the LVLH frame (rad/s)
    Omega_dot   - angular acceleration of the LVLH frame (rad/s^2)
    aA, aB      - absolute accelerations of A and B (km/s^2)
    r_rel       - position of B relative to A in GEF (km)
    v_rel       - velocity of B relative to A in GEF (km/s)
    a_rel       - acceleration of B relative to A in GEF (km/s^2)
    r_rel_x     - position of B relative to A in the LVLH frame
    v_rel_x     - velocity of B relative to A in the LVLH frame
    a_rel_x     - acceleration of B relative to A in the LVLH frame

    User h-functions required: None
*/
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> rva_relative(
    const std::vector<double>& rA,
    const std::vector<double>& vA,
    const std::vector<double>& rB,
    const std::vector<double>& vB,
    double mu)
{
    //...Calculate the vector hA:
    std::vector<double> hA = 
    {
        rA[1] * vA[2] - rA[2] * vA[1],
        rA[2] * vA[0] - rA[0] * vA[2],
        rA[0] * vA[1] - rA[1] * vA[0]
    };

    //...Calculate the unit vectors i, j and k:
    double rA_norm = std::sqrt(rA[0]*rA[0] + rA[1]*rA[1] + rA[2]*rA[2]);
    std::vector<double> i = {rA[0]/rA_norm, rA[1]/rA_norm, rA[2]/rA_norm};
    double hA_norm = std::sqrt(hA[0]*hA[0] + hA[1]*hA[1] + hA[2]*hA[2]);
    std::vector<double> k = {hA[0]/hA_norm, hA[1]/hA_norm, hA[2]/hA_norm};
    std::vector<double> j = 
    {
        k[1]*i[2] - k[2]*i[1],
        k[2]*i[0] - k[0]*i[2],
        k[0]*i[1] - k[1]*i[0]
    };

    //...Calculate the transformation matrix Qxx:
    std::vector<std::vector<double>> QXx = { i, j, k };

    //...Calculate Omega and Omega_dot:
    std::vector<double> Omega = 
    {
        hA[0] / (rA_norm * rA_norm),
        hA[1] / (rA_norm * rA_norm),
        hA[2] / (rA_norm * rA_norm)
    }; // Equation 7.5

    double dot_rA_vA = rA[0]*vA[0] + rA[1]*vA[1] + rA[2]*vA[2];
    double Omega_dot_scalar = -2.0 * (dot_rA_vA / (rA_norm * rA_norm));

    std::vector<double> Omega_dot = 
    {
        Omega_dot_scalar * Omega[0],
        Omega_dot_scalar * Omega[1],
        Omega_dot_scalar * Omega[2]
    }; // Equation 7.6

    //...Calculate the accelerations aA and aB:
    std::vector<double> aA = 
    {
        -mu * rA[0] / (rA_norm*rA_norm*rA_norm),
        -mu * rA[1] / (rA_norm*rA_norm*rA_norm),
        -mu * rA[2] / (rA_norm*rA_norm*rA_norm)
    };

    double rB_norm = std::sqrt(rB[0]*rB[0] + rB[1]*rB[1] + rB[2]*rB[2]);
    std::vector<double> aB = 
    {
        -mu * rB[0] / (rB_norm*rB_norm*rB_norm),
        -mu * rB[1] / (rB_norm*rB_norm*rB_norm),
        -mu * rB[2] / (rB_norm*rB_norm*rB_norm)
    };

    //...Calculate r_rel:
    std::vector<double> r_rel = 
    {
        rB[0] - rA[0],
        rB[1] - rA[1],
        rB[2] - rA[2]
    };

    //...Calculate v_rel:
    std::vector<double> v_rel = 
    {
        vB[0] - vA[0] - (Omega[1]*r_rel[2] - Omega[2]*r_rel[1]),
        vB[1] - vA[1] - (Omega[2]*r_rel[0] - Omega[0]*r_rel[2]),
        vB[2] - vA[2] - (Omega[0]*r_rel[1] - Omega[1]*r_rel[0])
    };

    //...Calculate a_rel:
    std::vector<double> cross_Omega_r_rel = {
        Omega[1]*r_rel[2] - Omega[2]*r_rel[1],
        Omega[2]*r_rel[0] - Omega[0]*r_rel[2],
        Omega[0]*r_rel[1] - Omega[1]*r_rel[0]
    };
    std::vector<double> cross_Omega_v_rel = {
        Omega[1]*v_rel[2] - Omega[2]*v_rel[1],
        Omega[2]*v_rel[0] - Omega[0]*v_rel[2],
        Omega[0]*v_rel[1] - Omega[1]*v_rel[0]
    };
    std::vector<double> cross_Omega_cross_Omega_r_rel = {
        Omega[1]*cross_Omega_r_rel[2] - Omega[2]*cross_Omega_r_rel[1],
        Omega[2]*cross_Omega_r_rel[0] - Omega[0]*cross_Omega_r_rel[2],
        Omega[0]*cross_Omega_r_rel[1] - Omega[1]*cross_Omega_r_rel[0]
    };
    std::vector<double> cross_Omega_dot_r_rel = {
        Omega_dot[1]*r_rel[2] - Omega_dot[2]*r_rel[1],
        Omega_dot[2]*r_rel[0] - Omega_dot[0]*r_rel[2],
        Omega_dot[0]*r_rel[1] - Omega_dot[1]*r_rel[0]
    };

    std::vector<double> a_rel = {
        aB[0] - aA[0]
          - cross_Omega_dot_r_rel[0]
          - cross_Omega_cross_Omega_r_rel[0]
          - 2.0 * cross_Omega_v_rel[0],
        aB[1] - aA[1]
          - cross_Omega_dot_r_rel[1]
          - cross_Omega_cross_Omega_r_rel[1]
          - 2.0 * cross_Omega_v_rel[1],
        aB[2] - aA[2]
          - cross_Omega_dot_r_rel[2]
          - cross_Omega_cross_Omega_r_rel[2]
          - 2.0 * cross_Omega_v_rel[2]
    };

    //...Calculate r_rel_x, v_rel_x and a_rel_x:
    std::vector<double> r_rel_x(3, 0.0), v_rel_x(3, 0.0), a_rel_x(3, 0.0);
    for (size_t iRow = 0; iRow < 3; ++iRow)
    {
        for (size_t jCol = 0; jCol < 3; ++jCol)
        {
            r_rel_x[iRow] += QXx[iRow][jCol] * r_rel[jCol];
            v_rel_x[iRow] += QXx[iRow][jCol] * v_rel[jCol];
            a_rel_x[iRow] += QXx[iRow][jCol] * a_rel[jCol];
        }
    }

    return {r_rel_x, v_rel_x, a_rel_x};
}

#endif // RVA_RELATIVE_H