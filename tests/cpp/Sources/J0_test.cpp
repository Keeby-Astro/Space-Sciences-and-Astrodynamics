#include <iostream>
#include <iomanip>
#include "J0.h"

/*
    This program computes J0 and the Julian day number using the data.
    
    year   - range: 1901 - 2099
    month  - range: 1 - 12
    day    - range: 1 - 31
    hour   - range: 0 - 23 (Universal Time)
    minute - rage: 0 - 60
    second - range: 0 - 60
    ut     - universal time (hr)
    j0     - Julian day number at 0 hr UT
    jd     - Julian day number at specified UT

    User h-function required: J0
*/
int main() 
{
    //...Data declaration:
    int year = 2004;
    int month = 5;
    int day = 12;

    int hour = 14;
    int minute = 45;
    int second = 30;
    //...

    double ut = hour + minute / 60.0 + second / 3600.0;

    //...Equation 5.46:
    double j0 = J0(year, month, day);

    //...Equation 5.47:
    double jd = j0 + ut / 24.0;

    //...Echo the input data and output the results to the command window:
    std::cout << "-----------------------------------------------------\n";
    std::cout << " Julian day calculation\n";
    std::cout << "\n Input data:\n";
    std::cout << "   Year            = " << year << "\n";
    std::cout << "   Month           = " << month << "\n";
    std::cout << "   Day             = " << day << "\n";
    std::cout << "   Hour            = " << hour << "\n";
    std::cout << "   Minute          = " << minute << "\n";
    std::cout << "   Second          = " << second << "\n";

    std::cout << "\n Julian day number = " << std::fixed << std::setprecision(3) << jd << "\n";
    std::cout << "-----------------------------------------------------\n";

    return 0;
}
