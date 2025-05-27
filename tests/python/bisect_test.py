import bisect

def bisect_test():
    '''
    This program uses the bisection method to find the three roots of
    Equation 2.204 for the earth-moon system.

    m1  - mass of the earth (kg)
    m2  - mass of the moon (kg)
    r12 - distance from the earth to the moon (km)
    p   - ratio of moon mass to total mass
    xl  - list containing the low-side estimates of the three roots
    xu  - list containing the high-side estimates of the three roots
    x   - list containing the three computed roots

    User py-function required: bisect
    User subfunction required: fun
    '''
    #...Input data:
    m1 = 5.974e24
    m2 = 7.348e22
    r12 = 3.844e5

    xl = [-1.1, 0.5, 1.0]
    xu = [-0.9, 1.0, 1.5]
    #...End input data

    p = m2 / (m1 + m2)

    x = []
    for i in range(3):
        x.append(bisect.bisect(lambda z: fun(z, p), xl[i], xu[i]))

    #...Output the results
    output(m1, m2, r12, p, x)

# -----------------
def fun(z, p):
# -----------------
    '''
    This subroutine evaluates the function in Equation 2.204

    z - the dimensionless x-coordinate
    p - defined above
    f - the value of the function
    '''
    f = (1 - p) * (z + p) / abs(z + p)**3 + p * (z + p - 1) / abs(z + p - 1)**3 - z
    return f

def output(m1, m2, r12, p, x):
    '''
    This function prints out the x-coordinates of L1, L2, and L3
    relative to the center of mass.
    '''
    #...Output to the console:
    print("\n\n---------------------------------------------")
    print("\n For:")
    print(f" m1  = {m1:.6g} kg")
    print(f" m2  = {m2:.6g} kg")
    print(f" r12 = {r12:.6g} km\n")
    print(" the 3 colinear Lagrange points (the roots of")
    print(" Equation 2.204) are:")
    print(f"\n L3: x = {x[0] * r12:10.6g} km (f(x3) = {fun(x[0], p):g})")
    print(f" L1: x = {x[1] * r12:10.6g} km (f(x1) = {fun(x[1], p):g})")
    print(f" L2: x = {x[2] * r12:10.6g} km (f(x2) = {fun(x[2], p):g})")
    print("\n---------------------------------------------\n")

if __name__ == "__main__":
    bisect_test()
