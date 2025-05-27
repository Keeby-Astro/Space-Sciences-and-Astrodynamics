// ALGORITHM 4.5: CALCULATION OF THE STATE VECTOR FROM THE ORBITAL ELEMENTS

#ifndef SV_FROM_COE_H
#define SV_FROM_COE_H

#include <cmath>
#include <vector>
#include <utility>

/*
    This function computes the state vector (r,v) from the
    classical orbital elements (coe).

    mu   - gravitational parameter (km^3/s^2)
    coe  - orbital elements [h e RA incl w TA]
           where
               h = angular momentum (km^2/s)
               e = eccentricity
               RA = right ascension of the ascending node (rad)
               incl = inclination of the orbit (rad)
               w = argument of perigee (rad)
               TA = true anomaly (rad)
    R3_w - Rotation matrix about the z-axis through the angle w
    R1_i - Rotation matrix about the x-axis through the angle i
    R3_W - Rotation matrix about the z-axis through the angle RA
    Q_pX - Matrix of the transformation from perifocal to geocentric
           equatorial frame
    rp   - position vector in the perifocal frame (km)
    vp   - velocity vector in the perifocal frame (km/s)
    r    - position vector in the geocentric equatorial frame (km)
    v    - velocity vector in the geocentric equatorial frame (km/s)
    
    User M-functions required: none
*/
inline void multiply_3x3(const double A[3][3], const double B[3][3], double C[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            C[i][j] = 0.0;
            for (int k = 0; k < 3; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

inline void transpose_3x3(const double A[3][3], double B[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            B[j][i] = A[i][j];
        }
    }
}


inline void multiply_3x3_by_3x1(const double A[3][3], const std::vector<double>& v, std::vector<double>& result)
{
    for (int i = 0; i < 3; i++)
    {
        result[i] = A[i][0] * v[0] + A[i][1] * v[1] + A[i][2] * v[2];
    }
}

inline std::pair<std::vector<double>, std::vector<double>> sv_from_coe(const std::vector<double>& coe, double mu)
{
    double h    = coe[0];
    double e    = coe[1];
    double RA   = coe[2];
    double incl = coe[3];
    double w    = coe[4];
    double TA   = coe[5];

    //...Equations 4.45 and 4.46 (rp and vp are column vectors):
    double factor = (h * h / mu) * (1.0 / (1.0 + e * std::cos(TA)));
    std::vector<double> rp(3);
    rp[0] = factor * std::cos(TA);
    rp[1] = factor * std::sin(TA);
    rp[2] = 0.0;

    double factor2 = (mu / h);
    std::vector<double> vp(3);
    vp[0] = factor2 * (-std::sin(TA));
    vp[1] = factor2 * ( e + std::cos(TA));
    vp[2] = 0.0;

    //...Equation 4.34:
    double R3_W[3][3] = {
        {  std::cos(RA),  std::sin(RA), 0.0 },
        { -std::sin(RA),  std::cos(RA), 0.0 },
        {           0.0,           0.0, 1.0 }
    };

    //...Equation 4.32:
    double R1_i[3][3] = {
        { 1.0,            0.0,            0.0 },
        { 0.0,  std::cos(incl),  std::sin(incl) },
        { 0.0, -std::sin(incl),  std::cos(incl) }
    };

    //...Equation 4.34:
    double R3_w[3][3] = {
        {  std::cos(w),  std::sin(w), 0.0 },
        { -std::sin(w),  std::cos(w), 0.0 },
        {          0.0,          0.0, 1.0 }
    };

    //...Equation 4.49:
    double M1[3][3], M2[3][3], Q_pX[3][3];
    multiply_3x3(R3_w, R1_i, M1);
    multiply_3x3(M1, R3_W, M2);
    transpose_3x3(M2, Q_pX);

    //...Equations 4.51 (r and v are column vectors):
    std::vector<double> r(3), v(3);
    multiply_3x3_by_3x1(Q_pX, rp, r);
    multiply_3x3_by_3x1(Q_pX, vp, v);

    //...Convert r and v into row vectors:
    return std::make_pair(r, v);
}

#endif  // SV_FROM_COE_H
