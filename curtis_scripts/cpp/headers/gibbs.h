// ALGORITHM 5.1: GIBBS' METHOD OF PRELIMINARY ORBIT DETERMINATION

#ifndef GIBBS_H
#define GIBBS_H

#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

/*
    This function uses the Gibbs method of orbit determination to
    to compute the velocity corresponding to the second of three
    supplied position vectors.

    mu              - gravitational parameter (km^3/s^2)
    R1, R2, R3      - three coplanar geocentric position vectors (km)
    r1, r2, r3      - the magnitudes of R1, R2 and R3 (km)
    c12, c23, c31   - three independent cross products among
                      R1, R2 and R3
    N, D, S         - vectors formed from R1, R2 and R3 during
                      the Gibbs  procedure
    tol             - tolerance for determining if R1, R2 and R3
                      are coplanar
    ierr            - = 0 if R1, R2, R3 are found to be coplanar
                      = 1 otherwise
    V2              - the velocity corresponding to R2 (km/s)

    User h-functions required: none
*/
const double mu = 398600.4418; // km^3/s^2

double norm(const std::array<double, 3>& vec) 
{
    return std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

std::array<double, 3> cross(const std::array<double, 3>& a, const std::array<double, 3>& b) 
{
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

double dot(const std::array<double, 3>& a, const std::array<double, 3>& b) 
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::array<double, 3> operator*(double scalar, const std::array<double, 3>& vec) 
{
    return {scalar * vec[0], scalar * vec[1], scalar * vec[2]};
}

std::array<double, 3> operator+(const std::array<double, 3>& a, const std::array<double, 3>& b) 
{
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

std::array<double, 3> operator*(const std::array<double, 3>& vec, double scalar) 
{
    return {vec[0] * scalar, vec[1] * scalar, vec[2] * scalar};
}

std::array<double, 3> operator/(const std::array<double, 3>& vec, double scalar) 
{
    return {vec[0] / scalar, vec[1] / scalar, vec[2] / scalar};
}

std::array<double, 3> operator-(const std::array<double, 3>& a, const std::array<double, 3>& b) 
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

std::pair<std::array<double, 3>, int> gibbs(const std::array<double, 3>& R1, const std::array<double, 3>& R2, const std::array<double, 3>& R3) {
    const double tol = 1e-4;
    int ierr = 0;

    //...Magnitudes of R1, R2 and R3:
    double r1 = norm(R1);
    double r2 = norm(R2);
    double r3 = norm(R3);

    //...Cross products among R1, R2 and R3:
    auto c12 = cross(R1, R2);
    auto c23 = cross(R2, R3);
    auto c31 = cross(R3, R1);

    //...Check that R1, R2 and R3 are coplanar; if not set error flag:
    if (std::abs(dot(R1, c23) / (r1 * norm(c23))) > tol) {
        ierr = 1;
        return {{0.0, 0.0, 0.0}, ierr};
    }

    //...Equation 5.13:
    auto N = r1 * c23 + r2 * c31 + r3 * c12;

    //...Equation 5.14:
    auto D = c12 + c23 + c31;

    //...Equation 5.21:
    auto S = R1 * (r2 - r3) + R2 * (r3 - r1) + R3 * (r1 - r2);

    //...Equation 5.22:
    auto V2 = std::sqrt(mu / norm(N) / norm(D)) * (cross(D, R2) / r2 + S);

    return {V2, ierr};
}

#endif // GIBBS_H
