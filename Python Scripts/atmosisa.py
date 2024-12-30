import numpy as np

def atmosisa(h_m):
    '''
    Returns standard atmosphere properties using an approximation
    to the International Standard Atmosphere (ISA) model for a
    given geometric altitude h_m (meters).

    T   - Temperature (K)
    a   - Speed of sound (m/s)
    P   - Pressure (Pa)
    rho - Density (g/m^3)

    User py-function required: none
    '''
    g0    = 9.80665
    R     = 287.053
    gamma = 1.4 

    # Layer definitions (U.S. Standard Atmosphere 1976 / ICAO up to ~85 km)
    layers = [
        {
            'h_base': 0.0,       'h_top': 11000.0,
            'T_base': 288.15,    # K
            'P_base': 101325.0,  # Pa
            'L':     -0.0065     # K/m
        },
        {
            'h_base': 11000.0,   'h_top': 20000.0,
            'T_base': 216.65,    # K
            'P_base': 22632.06,  # Pa
            'L':      0.0
        },
        {
            'h_base': 20000.0,   'h_top': 32000.0,
            'T_base': 216.65,    # K
            'P_base': 5474.89,   # Pa
            'L':      0.0010     # K/m
        },
        {
            'h_base': 32000.0,   'h_top': 47000.0,
            'T_base': 228.65,    # K
            'P_base': 868.02,    # Pa
            'L':      0.0028
        },
        {
            'h_base': 47000.0,   'h_top': 51000.0,
            'T_base': 270.65,    # K
            'P_base': 110.91,    # Pa
            'L':      0.0
        },
        {
            'h_base': 51000.0,   'h_top': 71000.0,
            'T_base': 270.65,    # K
            'P_base': 66.94,     # Pa
            'L':     -0.0028
        },
        {
            'h_base': 71000.0,   'h_top': 84852.0,
            'T_base': 214.65,    # K
            'P_base': 3.96,      # Pa
            'L':     -0.0020
        },
    ]

    # Clip altitudes to valid range
    if h_m < 0.0:
        raise ValueError("Altitude cannot be negative for this standard model.")
    if h_m > 84852.0:
        h_m = 84852.0

    # Find the layer in which the altitude resides
    layer_index = None
    for i, layer in enumerate(layers):
        if (i == len(layers) - 1) or (h_m <= layer['h_top']):
            layer_index = i
            break

    # Retrieve the layer data
    layer = layers[layer_index]
    h_b = layer['h_base']
    T_b = layer['T_base']
    P_b = layer['P_base']
    L   = layer['L']

    # Calculate temperature and pressure for altitude h_m in that layer
    delta_h = h_m - h_b
    if np.abs(L) > 1e-10:
        T = T_b + L * delta_h
        exponent = (-g0 / (R * L))
        P = P_b * (T / T_b) ** exponent
    else:
        T = T_b
        P = P_b * np.exp(-g0 * delta_h / (R * T_b))

    # Density and Speed of Sound
    rho = P / (R * T)
    a   = np.sqrt(gamma * R * T)

    print("Standard Atmosphere Properties")
    print("------------------------------") 
    print(f"Altitude: {h_m} m")
    print(f"Layer: {layer_index}")
    print(f"Temperature: {T} K")
    print(f"Pressure: {P} Pa")
    print(f"Density: {rho} kg/m^3")
    print(f"Speed of Sound: {a} m/s")
    print("------------------------------")

    return (T, a, P, rho)

if __name__ == "__main__":
    for h_m in [5000.0, 20000.0, 50000.0, 80000.0]:
        atmosisa(h_m)
