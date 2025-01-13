#include <iostream>
#include <iomanip>
#include "LST.h"

/*
    This program uses Algorithm 5.3 to obtain the local sidereal
    time from the data provided.

    lst   - local sidereal time (degrees)
    EL    - east longitude of the site (west longitude is negative):
               degrees (0 - 360)
               minutes (0 - 60)
               seconds (0 - 60)
    WL    - west longitude
    year  - range: 1901 - 2099
    month - range: 1 - 12
    day   - range: 1 - 31
    ut    - universal time
               hour (0 - 23)
               minute (0 - 60)
               second (0 - 60)

    User h-function required: LST
*/
int main() 
{
    //...Data declaration:
    //   East longitude:
    int degrees = 139;
    int minutes = 47;
    int seconds = 0;

    //   Date:
    int year = 2004;
    int month = 3;
    int day = 3;

    //   Universal time:
    int hour = 4;
    int minute = 30;
    int second = 0;
    //...

    //...Convert negative (west) longitude to east longitude:
    if (degrees < 0) 
    {
        degrees += 360;
    }

    //...Express the longitudes as decimal numbers:
    double EL = degrees + minutes / 60.0 + seconds / 3600.0;
    double WL = 360.0 - EL;

    //...Express universal time as a decimal number:
    double ut = hour + minute / 60.0 + second / 3600.0;

    //...Algorithm 5.3:
    double lst = LST(year, month, day, ut, EL);

    //......Echo the input data and output the results to the console:
    std::cout << "-----------------------------------------------------\n";
    std::cout << " Input data:\n";
    std::cout << "\n Year                      = " << year << "\n";
    std::cout << " Month                     = " << month << "\n";
    std::cout << " Day                       = " << day << "\n";
    std::cout << " UT (hr)                   = " << ut << "\n";
    std::cout << " West Longitude (deg)      = " << WL << "\n";
    std::cout << " East Longitude (deg)      = " << EL << "\n";
    std::cout << "\n Solution:\n";
    std::cout << "\n Local Sidereal Time (deg) = " << lst << "\n";
    std::cout << " Local Sidereal Time (hr)  = " << lst / 15.0 << "\n";
    std::cout << "-----------------------------------------------------\n";

    return 0;
}
