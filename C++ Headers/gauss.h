// ALGORITHMS 5.5 AND 5.6: GAUSS' METHOD OF PRELIMINARY ORBIT
// DETERMINATION WITH ITERATIVE IMPROVEMENT
#ifndef GAUSS_H
#define GAUSS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "kepler_U.h"
#include "f_and_g.h"
/*
    This function uses the Gauss method with iterative improvement
    (Algorithms 5.5 and 5.6) to calculate the state vector of an
    orbiting body from angles-only observations at three
    closely spaced times.

    mu                  - the gravitational parameter (km^3/s^2)
    t1, t2, t3          - the times of the observations (s)
    tau, tau1, tau3     - time interXs between observations (s)
    R1, R2, R3          - the observation site position vectors
                          at t1, t2, t3 (km)
    Rho1, Rho2, Rho3    - the direction cosine vectors of the
                          satellite at t1, t2, t3
    p1, p2, p3          - cross products among the three direction
                          cosine vectors
    Do                  - scalar triple product of Rho1, Rho2 and Rho3
    D                   - Matrix of the nine scalar triple products
                          of R1, R2 and R3 with p1, p2 and p3
    E                   - dot product of R2 and Rho2
    A, B                - constants in the expression relating slant range
                          to geocentric radius
    a,b,c               - coefficients of the 8th order polynomial
                          in the estimated geocentric radius x
    x                   - positive root of the 8th order polynomial
    rho1, rho2, rho3    - the slant ranges at t1, t2, t3
    r1, r2, r3          - the position vectors at t1, t2, t3 (km)
    r_old, v_old        - the estimated state vector at the end of
                          Algorithm 5.5 (km, km/s)
    rho1_old,
    rho2_old, and
    rho3_old            - the Xues of the slant ranges at t1, t2, t3
                          at the beginning of iterative improvement
                          (Algorithm 5.6) (km)
    diff1, diff2,
    and diff3           - the magnitudes of the differences between the
                          old and new slant ranges at the end of
                          each iteration
    tol                 - the error tolerance determining
                          convergence
    n                   - number of passes through the
                          iterative improvement loop
    nmax                - limit on the number of iterations
    ro, vo              - magnitude of the position and
                          velocity vectors (km, km/s)
    vro                 - radial velocity component (km)
    a                   - reciprocal of the semimajor axis (1/km)
    v2                  - computed velocity at time t2 (km/s)
    r, v                - the state vector at the end of Algorithm 5.6
                          (km, km/s)

    User h-functions required: kepler_U, f_and_g
    User subfunctions required: posroot
*/
inline std::vector<double> cross_product(const std::vector<double>& a,
                                         const std::vector<double>& b)
{
    return 
    {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

inline double dot_product(const std::vector<double>& a,
                          const std::vector<double>& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double posroot(const std::vector<double>& roots)
{
    /*
        This subfunction extracts the positive real roots from
        those obtained in the call to MATLAB's 'roots' function.
        If there is more than one positive root, the user is
        prompted to select the one to use.
    
        x        - the determined or selected positive root
        Roots    - the vector of roots of a polynomial
        posroots - vector of positive roots
    
        User h-functions required: none
    */
    std::vector<double> positive_roots;
    for (double X : roots)
    {
        if (X > 0.0)
            positive_roots.push_back(X);
    }

    if (positive_roots.empty())
    {
        throw std::runtime_error("No positive real roots found.");
    }
    return positive_roots.front();
}

inline std::vector<double> solve_polynomial_8th(const std::vector<double>& coeffs)
{
    double a = coeffs[2];
    double b = coeffs[5];
    double c = coeffs[8];

    auto f = [&](double x) {
        return std::pow(x, 8) + a * std::pow(x, 6) + b * std::pow(x, 3) + c;
    };

    double x_low = 1e-12;
    double x_high = 1e-2;

    const int MAX_EXPAND = 100;
    int expandCount = 0;

    double f_low = f(x_low);
    double f_high = f(x_high);

    while ((f_low * f_high > 0.0) && (expandCount < MAX_EXPAND))
    {
        x_low = x_high;
        x_high *= 10.0;
        f_low = f(x_low);
        f_high = f(x_high);
        expandCount++;
    }

    if (f_low == 0.0)
    {
        return {x_low};
    }
    if (f_high == 0.0)
    {
        return {x_high};
    }

    if (f_low * f_high > 0.0)
    {
        throw std::runtime_error(
            "Unable to bracket a positive real root for x^8 + a x^6 + b x^3 + c = 0."
        );
    }

    const int MAX_BISECT_ITER = 200;
    const double TOL = 1e-12;
    for (int i = 0; i < MAX_BISECT_ITER; ++i)
    {
        double x_mid = 0.5 * (x_low + x_high);
        double f_mid = f(x_mid);

        if (std::fabs(f_mid) < 1e-14)
        {
            return {x_mid};
        }

        if (f_low * f_mid > 0.0)
        {
            x_low = x_mid;
            f_low = f_mid;
        }
        else
        {
            x_high = x_mid;
            f_high = f_mid;
        }

        if (std::fabs(x_high - x_low) < TOL)
        {
            double x_final = 0.5 * (x_low + x_high);
            return {x_final};
        }
    }

    double x_guess = 0.5 * (x_low + x_high);
    return {x_guess};
}

inline void gauss(const std::vector<double>& Rho1,
                  const std::vector<double>& Rho2,
                  const std::vector<double>& Rho3,
                  const std::vector<double>& R1,
                  const std::vector<double>& R2,
                  const std::vector<double>& R3,
                  double t1,
                  double t2,
                  double t3,
                  std::vector<double>& r,
                  std::vector<double>& v,
                  std::vector<double>& r_old,
                  std::vector<double>& v_old,
                  double mu = 398600.0)
{
    //...Equation 5.98:
    double tau1 = t1 - t2;
    double tau3 = t3 - t2;

    //...Equation 5.101:
    double tau = tau3 - tau1;

    //...Independent cross products among the direction cosine vectors:
    std::vector<double> p1 = cross_product(Rho2, Rho3);
    std::vector<double> p2 = cross_product(Rho1, Rho3);
    std::vector<double> p3 = cross_product(Rho1, Rho2);

    //...Equation 5.108:
    double Do = dot_product(Rho1, p1);

    //...Equations 5.109b, 5.110b, 5.111b:
    std::vector<std::vector<double>> D(3, std::vector<double>(3));
    D[0][0] = dot_product(R1, p1); D[0][1] = dot_product(R1, p2); D[0][2] = dot_product(R1, p3);
    D[1][0] = dot_product(R2, p1); D[1][1] = dot_product(R2, p2); D[1][2] = dot_product(R2, p3);
    D[2][0] = dot_product(R3, p1); D[2][1] = dot_product(R3, p2); D[2][2] = dot_product(R3, p3);

    //...Equation 5.115b:
    double E = dot_product(R2, Rho2);

    //...Equations 5.112b and 5.112c:
    double A = (1.0 / Do) * ( -D[0][1]*tau3/tau + D[1][1] + D[2][1]*tau1/tau );
    double B = (1.0 / (6.0 * Do)) *
               ( D[0][1] * ((tau3*tau3 - tau*tau) * (tau3/tau)) +
                 D[2][1] * ((tau*tau - tau1*tau1) * (tau1/tau)) );

    //...Equations 5.117:
    double R2norm = std::sqrt(dot_product(R2, R2));
    double a_coeff = -(A*A + 2.0*A*E + R2norm*R2norm);
    double b_coeff = -2.0*mu*B*(A + E);
    double c_coeff = -(mu*B)*(mu*B);

    //...Calculate the roots of Equation 5.116:
    std::vector<double> poly_coeffs = {1, 0, a_coeff, 0, 0, b_coeff, 0, 0, c_coeff};
    std::vector<double> poly_roots;
    try
    {
        poly_roots = solve_polynomial_8th(poly_coeffs);
    }
    catch (const std::runtime_error& e)
    {
        throw;
    }

    //...Find the positive real root:
    double x = posroot(poly_roots);

    //...Equations 5.99a and 5.99b:
    double f1 = 1.0 - 0.5 * mu * (tau1*tau1) / std::pow(x, 3);
    double f3 = 1.0 - 0.5 * mu * (tau3*tau3) / std::pow(x, 3);

    //...Equations 5.100a and 5.100b:
    double g1 = tau1 - (1.0/6.0)*mu*std::pow((tau1/x), 3);
    double g3 = tau3 - (1.0/6.0)*mu*std::pow((tau3/x), 3);

    //...Equation 5.112a:
    double rho2 = A + mu*B/std::pow(x, 3);

    //...Equations 5.113:
    double numerator_rho1 = 6.0 * ((D[2][0]*(tau1/tau3) + D[1][0]*(tau/tau3))*std::pow(x,3))
                           + mu * D[2][0]*( (tau*tau) - (tau1*tau1) )*tau1/tau3;
    double denominator_rho1 = 6.0*std::pow(x,3) + mu*( (tau*tau) - (tau3*tau3) );
    double rho1 = (1.0 / Do) * ( ( numerator_rho1 / denominator_rho1 ) - D[0][0]);

    //...Equations 5.114:
    double numerator_rho3 = 6.0 * ((D[0][2]*(tau3/tau1) - D[1][2]*(tau/tau1))*std::pow(x,3))
                           + mu * D[0][2]*( (tau*tau) - (tau3*tau3) )*tau3/tau1;
    double denominator_rho3 = 6.0*std::pow(x,3) + mu*( (tau*tau) - (tau1*tau1) );
    double rho3 = (1.0 / Do) * ( ( numerator_rho3 / denominator_rho3 ) - D[2][2]);

    //...Equations 5.86:
    std::vector<double> r1(3), r2v(3), r3v(3);
    for (int i = 0; i < 3; ++i)
    {
        r1[i]  = R1[i] + rho1*Rho1[i];
        r2v[i] = R2[i] + rho2*Rho2[i];
        r3v[i] = R3[i] + rho3*Rho3[i];
    }

    //...Equation 5.118:
    std::vector<double> v2(3);
    double denom_v2 = (f1*g3 - f3*g1);
    for (int i = 0; i < 3; ++i)
    {
        v2[i] = ( -f3*r1[i] + f1*r3v[i] ) / denom_v2;
    }

    //...Save the initial estimates of r2 and v2:
    r_old = r2v;
    v_old = v2;

    //...End of Algorithm 5.5

    //...Use Algorithm 5.6 to improve the accuracy of the initial estimates.
    //...Initialize the iterative improvement loop and set error tolerance:
    double rho1_old = rho1;
    double rho2_old = rho2;
    double rho3_old = rho3;
    double diff1    = 1.0;
    double diff2    = 1.0;
    double diff3    = 1.0;
    int n    = 0;
    int nmax = 1000;
    double tol = 1.e-8;

    //...Iterative improvement loop:
    while ( (diff1 > tol) && (diff2 > tol) && (diff3 > tol) && (n < nmax) )
    {
        n++;

        //...Compute quantities required by universal Kepler's equation:
        double ro = std::sqrt(dot_product(r2v, r2v));
        double vo = std::sqrt(dot_product(v2, v2));
        double vro = dot_product(v2, r2v)/ro;
        double a_recip = 2.0/ro - (vo*vo)/mu;

        //...Solve universal Kepler's equation at times tau1 and tau3 for
        // universal anomalies x1 and x3:
        double x1 = kepler_U(tau1, ro, vro, a_recip, mu);
        double x3 = kepler_U(tau3, ro, vro, a_recip, mu);

        //...Calculate the Lagrange f and g coefficients at times tau1
        // and tau3:
        double ff1, gg1, ff3, gg3;
        f_and_g(x1, tau1, ro, a_recip, ff1, gg1, mu);
        f_and_g(x3, tau3, ro, a_recip, ff3, gg3, mu);

        //...Update the f and g functions at times tau1 and tau3 by
        // averaging old and new:
        f1 = (f1 + ff1)/2.0;
        f3 = (f3 + ff3)/2.0;
        g1 = (g1 + gg1)/2.0;
        g3 = (g3 + gg3)/2.0;

        //...Equations 5.96 and 5.97:
        double c1 = g3/(f1*g3 - f3*g1);
        double c3 = -g1/(f1*g3 - f3*g1);

        //...Equations 5.109a, 5.110a and 5.111a:
        rho1 = (1.0/Do)*( -D[0][0] + (1.0/c1)*D[1][0] - (c3/c1)*D[2][0] );
        rho2 = (1.0/Do)*( -c1*D[0][1] +       D[1][1] -  c3*D[2][1] );
        rho3 = (1.0/Do)*( -(c1/c3)*D[0][2] + (1.0/c3)*D[1][2] - D[2][2] );

        //...Equations 5.86:
        for (int i = 0; i < 3; ++i)
        {
            r1[i]  = R1[i] + rho1*Rho1[i];
            r2v[i] = R2[i] + rho2*Rho2[i];
            r3v[i] = R3[i] + rho3*Rho3[i];
        }

        //...Equation 5.118:
        denom_v2 = (f1*g3 - f3*g1);
        for (int i = 0; i < 3; ++i)
        {
            v2[i] = ( -f3*r1[i] + f1*r3v[i] ) / denom_v2;
        }

        diff1 = std::fabs(rho1 - rho1_old);
        diff2 = std::fabs(rho2 - rho2_old);
        diff3 = std::fabs(rho3 - rho3_old);

        //...Update the slant ranges:
        rho1_old = rho1;
        rho2_old = rho2;
        rho3_old = rho3;
    }

    std::cout << "\n( **Number of Gauss improvement iterations = " << n << " )\n\n";
    if (n >= nmax)
    {
        std::cout << "\n\n **Number of iterations exceeds " << nmax << " \n\n";
    }
    //...End iterative improvement loop

    // Return the state vector for the central observation:
    r = r2v;
    v = v2;
}

#endif // GAUSS_H
