#include <iostream>
#include <vector>
#include "ra_and_dec_from_r.h"

/*
    This program calculates the right ascension and declination
    from the geocentric equatorial position vector using the data.

    r   - position vector r (km)
    ra  - right ascension (deg)
    dec - declination (deg)

    User h-functions required: ra_and_dec_from_r
*/

int main()
{
    std::vector<double> r = {-5368.0, -1784.0, 3691.0};
    auto [ra, dec] = ra_and_dec_from_r(r);

    std::cout << "\n -----------------------------------------------------\n";
    std::cout << "\n r               = [" << r[0] << " " << r[1] << " " << r[2] << "] (km)";
    std::cout << "\n right ascension = " << ra << " deg";
    std::cout << "\n declination     = " << dec << " deg";
    std::cout << "\n\n -----------------------------------------------------\n";

    return 0;
}
