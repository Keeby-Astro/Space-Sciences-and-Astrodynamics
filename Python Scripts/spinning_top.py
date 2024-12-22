import numpy as np
import matplotlib.pyplot as plt
import rkf45
import q_from_dcm
import dcm_from_q
import dcm_to_euler

def spinning_top():
    '''
    This program numerically integrates Euler's equations of motion
    for the spinning top. The  quaternion is used to obtain the time history
    of the top's orientation.
    
    User py-functions required: rkf45, q_from_dcm, dcm_from_q, dcm_to_euler
    User subfunction required: rates
    '''
    # Data declaration:
    g = 9.807  # Acceleration of gravity (m/s^2)
    m = 0.5  # Mass in kg
    d = 0.05  # Distance of center of mass from pivot point (m)
    A = 12.e-4  # Moment of inertia about body x (kg-m^2)
    B = 12.e-4  # Moment of inertia about body y (kg-m^2)
    C = 4.5e-4  # Moment of inertia about body z (kg-m^2)
    ws0 = 1000 * 2 * np.pi / 60  # Spin rate (rad/s)

    wp0 = 51.93 * 2 * np.pi / 60  # Precession rate (rad/s) Use to obtain Fig. 11.33
    wp0 = 0  # Use to obtain Fig, 11.34

    wn0 = 0  # Nutation rate (deg/s)
    theta = 60  # Initial nutation angle (deg)
    z = np.array([0, -np.sin(np.radians(theta)), np.cos(np.radians(theta))])  # Initial z-axis direction:
    p = np.array([1, 0, 0])  # Initial x-axis direction (or a line defining x-z plane)

    y = np.cross(z, p)  # y-axis direction (normal to x-z plane)
    x = np.cross(y, z)  # x-axis direction (normal to y-z plane)
    i = x / np.linalg.norm(x)  # Unit vector along x axis
    j = y / np.linalg.norm(y)  # Unit vector along y axis
    k = z / np.linalg.norm(z)  # Unit vector along z axis
    QXx = np.array([i, j, k])  # Initial direction cosine matrix

    # Initial precession, nutation, and spin angles (deg):
    phi0, theta0, psi0 = dcm_to_euler.dcm_to_euler(QXx)

    # Initial quaternion (column vector):
    q0 = q_from_dcm.q_from_dcm(QXx)

    # Initial body-frame angular velocity, column vector (rad/s):
    w0 = np.array([
        wp0 * np.sin(np.radians(theta0)) * np.sin(np.radians(psi0)) + wn0 * np.cos(np.radians(psi0)),
        wp0 * np.sin(np.radians(theta0)) * np.cos(np.radians(psi0)) - wn0 * np.sin(np.radians(psi0)),
        ws0 + wp0 * np.cos(np.radians(theta0))
    ])

    t0 = 0  # Initial time (s)
    tf = 1.153  # Final time (s)  (for 360 degrees of precession)
    f0 = np.concatenate((q0, w0))  # Initial conditions vector (quaternion & angular velocities)

    # RKF4(5) numerical ODE solver. Time derivatives computed in function 'rates' below.
    t, f = rkf45.rkf45(rates, [t0, tf], f0)

    # Solutions for quaternion and angular velocities at 'nsteps' times from t0 to tf
    q = f[:, :4]
    wx = f[:, 4]
    wy = f[:, 5]
    wz = f[:, 6]

    # Obtain the direction cosine matrix, the Euler angles and the Euler angle rates at each solution time:
    prec, nut, spin, wp, wn, ws = [], [], [], [], [], []
    for m in range(len(t)):
        # DCM from the quaternion:
        QXx = dcm_from_q.dcm_from_q(q[m, :])
        # Euler angles (deg) from DCM:
        phi, theta, psi = dcm_to_euler.dcm_to_euler(QXx)
        prec.append(phi)
        nut.append(theta)
        spin.append(psi)
        # Euler rates from Eqs. 11.116:
        wp.append((wx[m] * np.sin(np.radians(psi)) + wy[m] * np.cos(np.radians(psi))) / np.sin(np.radians(theta)))
        wn.append(wx[m] * np.cos(np.radians(psi)) - wy[m] * np.sin(np.radians(psi)))
        ws.append(-wp[-1] * np.cos(np.radians(theta)) + wz[m])

    plotit(t, prec, wp, nut, wn, spin, ws)

def rates(t, f):
    g = 9.807  # Acceleration of gravity (m/s^2)
    m = 0.5  # Mass in kg
    d = 0.05  # Distance of center of mass from pivot point (m)
    A = 12.e-4  # Moment of inertia about body x (kg-m^2)
    B = 12.e-4  # Moment of inertia about body y (kg-m^2)
    C = 4.5e-4  # Moment of inertia about body z (kg-m^2)

    q = f[:4]  # components of quaternion
    wx, wy, wz = f[4:]  # angular velocity along x, y, z

    q = q / np.linalg.norm(q)  # normalize the quaternion

    Q = dcm_from_q.dcm_from_q(q)  # DCM from quaternion

    # Body frame components of the moment of the weight vector about the pivot point:
    M = Q @ np.array([-m * g * d * Q[2, 1], m * g * d * Q[2, 0], 0])

    # Skew-symmetric matrix of angular velocities:
    Omega = np.array([
        [0, wz, -wy, wx],
        [-wz, 0, wx, wy],
        [wy, -wx, 0, wz],
        [-wx, -wy, -wz, 0]
    ])
    q_dot = Omega @ q / 2  # time derivative of quaternion

    # Euler's equations:
    wx_dot = M[0] / A - (C - B) * wy * wz / A  # time derivative of wx
    wy_dot = M[1] / B - (A - C) * wz * wx / B  # time derivative of wy
    wz_dot = M[2] / C - (B - A) * wx * wy / C  # time derivative of wz

    # Return the rates in a column vector:
    return np.concatenate((q_dot, [wx_dot, wy_dot, wz_dot]))

def plotit(t, prec, wp, nut, wn, spin, ws):
    plt.figure('Euler angles and their rates', figsize=(10, 8))

    plt.subplot(321)
    plt.plot(t, prec)
    plt.xlabel('time (s)')
    plt.ylabel('Precession angle (deg)')
    plt.grid()

    plt.subplot(322)
    plt.plot(t, np.array(wp) * 60 / (2 * np.pi))
    plt.xlabel('time (s)')
    plt.ylabel('Precession rate (rpm)')
    plt.grid()

    plt.subplot(323)
    plt.plot(t, nut)
    plt.xlabel('time (s)')
    plt.ylabel('Nutation angle (deg)')
    plt.grid()

    plt.subplot(324)
    plt.plot(t, np.array(wn) * 180 / np.pi)
    plt.xlabel('time (s)')
    plt.ylabel('Nutation rate (deg/s)')
    plt.grid()

    plt.subplot(325)
    plt.plot(t, spin)
    plt.xlabel('time (s)')
    plt.ylabel('Spin angle (deg)')
    plt.grid()

    plt.subplot(326)
    plt.plot(t, np.array(ws) * 60 / (2 * np.pi))
    plt.xlabel('time (s)')
    plt.ylabel('Spin rate (rpm)')
    plt.grid()

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    spinning_top()