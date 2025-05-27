// ALGORITHM 2.4: FIND THE ROOT OF A FUNCTION USING THE
// BISECTION METHOD

#ifndef BISECT_H
#define BISECT_H

#include <cmath>
#include <functional>

/*
    This function evaluates a root of a function using
    the bisection method.

    tol  - error to within which the root is computed
    n    - number of iterations
    xl   - low end of the interval containing the root
    xu   - upper end of the interval containing the root
    i    - loop index
    xm   - mid-point of the interval from xl to xu
    fun  - name of the function whose root is being found
    fxl  - value of fun at xl
    fxm  - value of fun at xm
    root - the computed root

    User h-functions required: none
*/
double bisect(std::function<double(double)> fun, double xl, double xu, double tol = 1.e-6)
{
    int n = static_cast<int>(std::ceil(std::log(std::abs(xu - xl) / tol) / std::log(2)));
    double xm = 0.0;

    for (int i = 0; i < n; ++i)
    {
        xm = (xl + xu) / 2.0;
        double fxl = fun(xl);
        double fxm = fun(xm);
        if (fxl * fxm > 0)
        {
            xl = xm;
        }
        else
        {
            xu = xm;
        }
    }

    return xm;
}

#endif // BISECT_H