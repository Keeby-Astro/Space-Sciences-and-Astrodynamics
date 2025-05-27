import ra_and_dec_from_r

def ra_and_dec_from_r_test():
    '''
    This program calculates the right ascension and declination
    from the geocentric equatorial position vector using the data.

    r   - position vector r (km)
    ra  - right ascension (deg)
    dec - declination (deg)

    User py-functions required: ra_and_dec_from_r
    '''
    r = [-5368, -1784, 3691]
    ra, dec = ra_and_dec_from_r.ra_and_dec_from_r(r)

    print("\n -----------------------------------------------------")
    print(f"\n r               = [{r[0]} {r[1]} {r[2]}] (km)")
    print(f" right ascension = {ra} deg")
    print(f" declination     = {dec} deg")
    print("\n -----------------------------------------------------\n")

if __name__ == "__main__":
    ra_and_dec_from_r_test()