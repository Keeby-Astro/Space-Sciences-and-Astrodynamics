import J0

def J0_test():
    '''
    This program computes J0 and the Julian day number using the data.

    year   - range: 1901 - 2099
    month  - range: 1 - 12
    day    - range: 1 - 31
    hour   - range: 0 - 23 (Universal Time)
    minute - range: 0 - 60
    second - range: 0 - 60
    ut     - universal time (hr)
    j0     - Julian day number at 0 hr UT
    jd     - Julian day number at specified UT

    User py-function required: J0
    '''
    #...Data declaration:
    year = 2004
    month = 5
    day = 12

    hour = 14
    minute = 45
    second = 30
    #...

    ut = hour + minute / 60 + second / 3600

    #...Equation 5.46:
    j0 = J0.J0(year, month, day)

    #...Equation 5.47:
    jd = j0 + ut / 24

    #...Echo the input data and output the results to the console:
    print('-----------------------------------------------------')
    print('\n Example 5.4: Julian day calculation\n')
    print('\n Input data:\n')
    print(f'   Year            = {year}')
    print(f'   Month           = {month}')
    print(f'   Day             = {day}')
    print(f'   Hour            = {hour}')
    print(f'   Minute          = {minute}')
    print(f'   Second          = {second}\n')

    print(f' Julian day number = {jd:11.3f}')
    print('-----------------------------------------------------')

if __name__ == '__main__':
    J0_test()