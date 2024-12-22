import f_and_g

def f_and_g_test():
    '''
    This function tests the f_and_g function by providing sample inputs
    and verifying the outputs. 

    User py-functions required: f_and_g, stumpC, stumpS
    '''
    mu = 398600.4418 # km^3/s^2, standard for Earth

    #...Inputs
    ro = 7000  # km, radial position at time to
    t = 500  # s, time elapsed
    x = 1.5  # universal anomaly
    a = 1 / 10000  # 1/km, reciprocal of the semimajor axis

    # Call the f_and_g function
    f, g = f_and_g.f_and_g(x, t, ro, a, mu)

    #...Results
    print('Results')
    print(f'Computed f: {f:.6f}')
    print(f'Computed g: {g:.6f}')

if __name__ == "__main__":
    f_and_g_test()