// ALGORITHM 1.1: NUMERICAL INTEGRATION BY RUNGE-KUTTA METHODS RK1, RK2, RK3, OR RK4

#ifndef RK1_4_H
#define RK1_4_H

#include <iostream>
#include <vector>
#include <functional>
#include <stdexcept>

struct RKResult
{
    std::vector<double> tout;
    std::vector<std::vector<double>> yout;
};

/*
    This function uses a selected Runge-Kutta procedure to integrate
    a system of first-order differential equations dy/dt = f(t,y).
    
    y               - column vector of solutions
    f               - column vector of the derivatives dy/dt
    t               - time
    rk              - = 1 for RK1; = 2 for RK2; = 3 for RK3; = 4 for RK4
    n_stages        - the number of points within a time interval that
                      the derivatives are to be computed
    a               - coefficients for locating the solution points within
                      each time interval
    b               - coefficients for computing the derivatives at each
    interior point
    c               - coefficients for the computing solution at the end of
                      the time step
    ode_function    - handle for user M-function in which the derivatives f
                      are computed
    tspan           - the vector [t0 tf] giving the time interval for the
                      solution
    t0              - initial time
    tf              - final time
    y0              - column vector of initial values of the vector y
    tout            - column vector of times at which y was evaluated
    yout            - a matrix, each row of which contains the components of y
                      evaluated at the correponding time in tout
    h               - time step
    ti              - time at the beginning of a time step
    yi              - values of y at the beginning of a time step
    t_inner         - time within a given time step
    y_inner         - values of y within a given time step

    User M-function required: ode_function
*/

class RK1_4 
{
public:
    static RKResult integrate(
        std::function<std::vector<double>(double, const std::vector<double>&)> ode_function,
        const std::pair<double, double>& tspan,
        const std::vector<double>& y0,
        double h,
        int rk
    ) {
        int n_stages;
        std::vector<double> a, c;
        std::vector<std::vector<double>> b;

        //...Determine which of the four Runge-Kutta methods is to be used:
        switch (rk) {
            case 1:
                n_stages = 1;
                a = {0};
                b = {{0}};
                c = {1};
                break;
            case 2:
                n_stages = 2;
                a = {0, 1};
                b = {{0, 0}, {1, 0}};
                c = {0.5, 0.5};
                break;
            case 3:
                n_stages = 3;
                a = {0, 0.5, 1};
                b = {{0, 0, 0}, {0.5, 0, 0}, {-1, 2, 0}};
                c = {1.0 / 6, 2.0 / 3, 1.0 / 6};
                break;
            case 4:
                n_stages = 4;
                a = {0, 0.5, 0.5, 1};
                b = {{0, 0, 0, 0}, {0.5, 0, 0, 0}, {0, 0.5, 0, 0}, {0, 0, 1, 0}};
                c = {1.0 / 6, 1.0 / 3, 1.0 / 3, 1.0 / 6};
                break;
            default:
                throw std::invalid_argument("The parameter rk must have the value 1, 2, 3, or 4.");
        }

        double t0 = tspan.first;
        double tf = tspan.second;
        double t = t0;
        std::vector<double> y = y0;
        RKResult result;
        result.tout.push_back(t);
        result.yout.push_back(y);

        while (t < tf) {
            double ti = t;
            std::vector<double> yi = y;

            std::vector<std::vector<double>> f(n_stages, std::vector<double>(y0.size(), 0));

            //...Evaluate the time derivative(s) at the 'n_stages' points within the
            // current interval:
            for (int i = 0; i < n_stages; ++i) {
                double t_inner = ti + a[i] * h;
                std::vector<double> y_inner = yi;

                for (int j = 0; j < i; ++j) {
                    for (size_t k = 0; k < y0.size(); ++k) {
                        y_inner[k] += h * b[i][j] * f[j][k];
                    }
                }

                f[i] = ode_function(t_inner, y_inner);
            }

            //...Update the time step:
            h = std::min(h, tf - t);
            t += h;
            for (size_t k = 0; k < y.size(); ++k) {
                for (int i = 0; i < n_stages; ++i) {
                    y[k] += h * f[i][k] * c[i];
                }
            }

            result.tout.push_back(t); // adds t to the bottom of the column vector tout
            result.yout.push_back(y); // adds y to the bottom of the matrix yout
        }

        return result;
    }
};

#endif
