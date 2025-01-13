/*
ALGORITHM 1.3: NUMERICAL INTEGRATION OF A SYSTEM OF
FIRST-ORDER DIFFERENTIAL EQUATIONS BY THE RUNGE-KUTTA-FEHLBERG
4(5) METHOD WITH ADAPTIVE SIZE CONTROL
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <limits>
#include <algorithm>

/*
{
    This function uses the Runge-Kutta-Fehlberg 4(5) algorithm to
    integrate a system of first-order differential equations
    dy/dt = f(t,y).

    y               - column vector of solutions
    f               - column vector of the derivatives dy/dt
    t               - time
    a               - Fehlberg coefficients for locating the six solution
                      points (nodes) within each time interval.
    b               - Fehlberg coupling coefficients for computing the
                      derivatives at each interior point
    c4              - Fehlberg coefficients for the fourth-order solution
    c5              - Fehlberg coefficients for the fifth-order solution
    tol             - allowable truncation error
    ode_function    - handle (in C++, a std::function) for the user function
                      that computes the derivatives f
    tspan           - the pair [t0, tf] giving the time interval for the
                      solution
    t0              - initial time
    tf              - final time
    y0              - column vector of initial values of the vector y
    tout            - column vector of times at which y was evaluated
    yout            - a matrix, each row of which contains the components of y
                      evaluated at the corresponding time in tout
    h               - time step
    hmin            - minimum allowable time step
    ti              - time at the beginning of a time step
    yi              - values of y at the beginning of a time step
    t_inner         - time within a given time step
    y_inner         - values of y within a given time step
    te              - truncation error for each y at a given time step
    te_allowed      - allowable truncation error
    te_max          - maximum absolute value of the components of te
    ymax            - maximum absolute value of the components of y
    tol             - relative tolerance
    delta           - fractional change in step size
    eps             - unit roundoff error (the smallest number for which
                      1 + eps > 1)
    eps(x)          - the smallest number such that x + eps(x) = x
}
*/

struct RKF45 {
    std::vector<double>              tout;  // Times at which solutions are stored
    std::vector<std::vector<double>> yout;  // Solutions at each time
};

/*
-------------------------------------------------------------------------
*/

RKF45 rkf45(
    std::function<std::vector<double>(double, const std::vector<double>&)> ode_function,
    const std::pair<double,double>&   tspan,
    const std::vector<double>&        y0,
    double                            tolerance = 1.e-8)
{
    // Convert tspan to t0 and tf
    double t0 = tspan.first;
    double tf = tspan.second;

    // Initialize time and solution vector
    double t = t0;
    std::vector<double> y = y0;

    // Prepare output containers
    RKF45 result;
    result.tout.push_back(t);
    result.yout.push_back(y);

    // a, b, c4, c5 are the Fehlberg coefficients
    static const double a[6] = {0.0, 0.25, 3.0/8.0, 12.0/13.0, 1.0, 0.5};

    static const double b[6][5] = {
        {  0.0,       0.0,       0.0,       0.0,       0.0     },
        {  0.25,      0.0,       0.0,       0.0,       0.0     },
        {  3.0/32.0,  9.0/32.0,  0.0,       0.0,       0.0     },
        { 1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0,  0.0,  0.0 },
        {  439.0/216.0,   -8.0,  3680.0/513.0, -845.0/4104.0, 0.0 },
        {  -8.0/27.0,     2.0,  -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0}
    };

    static const double c4[6] = {25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0 };
    static const double c5[6] = {16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0 };

    // Store user tolerance
    double tol = tolerance;

    // Assumed initial step
    double h = (tf - t0) / 100.0;

    // The number of ODEs in the system
    std::size_t n = y0.size();

    std::vector<std::vector<double>> f(n, std::vector<double>(6, 0.0));

    // Main integration loop
    while (t < tf)
    {
        double eps_t = std::nextafter(t, t + 1.0) - t;
        double hmin  = 16.0 * eps_t;

        double ti = t;
        std::vector<double> yi = y;  // Copy current y

        // Evaluate the time derivatives at six points within the current interval
        for (int i = 0; i < 6; ++i)
        {
            double t_inner = ti + a[i] * h;
            std::vector<double> y_inner = yi;

            // y_inner = yi + h * sum_{j=0..i-1}[ b[i][j]*f(:,j) ]
            for (int j = 0; j < i; ++j)
            {
                for (std::size_t eq = 0; eq < n; ++eq)
                {
                    y_inner[eq] += h * b[i][j] * f[eq][j];
                }
            }

            // Compute derivative f(t_inner, y_inner)
            std::vector<double> f_temp = ode_function(t_inner, y_inner);
            for (std::size_t eq = 0; eq < n; ++eq)
            {
                f[eq][i] = f_temp[eq];
            }
        }

        // Compute the maximum truncation error
        std::vector<double> teVec(n, 0.0);
        for (std::size_t eq = 0; eq < n; ++eq)
        {
            double sum_c4_minus_c5 = 0.0;
            for (int k = 0; k < 6; ++k)
            {
                sum_c4_minus_c5 += f[eq][k] * (c4[k] - c5[k]);
            }
            teVec[eq] = h * sum_c4_minus_c5;
        }

        double te_max = 0.0;
        for (std::size_t eq = 0; eq < n; ++eq)
        {
            te_max = std::max(te_max, std::fabs(teVec[eq]));
        }

        // Compute the allowable truncation error
        double ymax = 0.0;
        for (std::size_t eq = 0; eq < n; ++eq)
        {
            ymax = std::max(ymax, std::fabs(y[eq]));
        }
        double te_allowed = tol * std::max(ymax, 1.0);

        // Compute the fractional change in step size
        double denom   = te_max + std::numeric_limits<double>::epsilon();
        double delta   = std::pow((te_allowed / denom), 0.2);

        // If the truncation error is in bounds, then update the solution
        if (te_max <= te_allowed)
        {
            // Prevent overshooting tf
            if (t + h > tf) {
                h = tf - t;
            }

            t = t + h;
            std::vector<double> newY(n, 0.0);
            for (std::size_t eq = 0; eq < n; ++eq)
            {
                double sum_c5 = 0.0;
                for (int k = 0; k < 6; ++k)
                {
                    sum_c5 += f[eq][k] * c5[k];
                }
                newY[eq] = yi[eq] + h * sum_c5;
            }
            y = newY;

            // Store output
            result.tout.push_back(t);
            result.yout.push_back(y);
        }

        // Update the time step
        h = std::min(delta * h, 4.0 * h);

        // If h < hmin, then step size fell below minimum allowable value
        if (h < hmin)
        {
            std::cerr << "\n\n Warning: Step size fell below its minimum"
                      << " allowable value (" << hmin << ") at time " << t << ".\n\n";
            break;
        }
    }

    return result;
}