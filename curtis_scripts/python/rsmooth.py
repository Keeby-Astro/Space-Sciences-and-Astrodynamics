import numpy as np
from scipy.fft import dct, idct
from scipy.optimize import minimize_scalar

def rsmooth(y):
    '''
    RSmooth - Robust smoothing of gridded data in one and higher dimensions.

    This function implements the algorithm described by:
    Garcia D. Robust smoothing of gridded data in one and higher dimensions with missing values.
    Comput Stat Data Anal. 2010 Apr 1;54(4):1167-1178. doi: 10.1016/j.csda.2009.09.020.
    PMID: 24795488; PMCID: PMC4008475.

    y - Input data array to be smoothed.
    z - Smoothed output array.

    Note: Please ensure the proper acknowledgment of the original author
            and source when using this code.
    '''
    y = np.asarray(y, dtype=float)
    if y.ndim != 1 or not np.isfinite(y).all():
        raise ValueError("Input must be a finite one-dimensional series.")
    n = len(y)
    if n < 3 or np.ptp(y) <= 100 * np.finfo(float).eps * max(1.0, np.max(np.abs(y))):
        return y.copy()

    center, scale = np.mean(y), np.ptp(y)
    data = (y - center) / scale
    Lambda = -2 + 2 * np.cos(np.arange(n) * np.pi / n)
    Lambda_sq = Lambda**2
    W = np.ones(n)
    z = data.copy()
    for k in range(6):
        for iteration in range(2000):
            DCTy = dct(W * (data - z) + z, norm='ortho')

            def GCVscore(p):
                Gamma = 1 / (1 + 10**p * Lambda_sq)
                fitted = idct(Gamma * DCTy, norm='ortho')
                RSS = np.sum(W * (data - fitted)**2)
                return RSS / n / max(1 - np.sum(Gamma) / n, 1.e-15)**2
            grid = np.linspace(-6, 16, 45)
            scores = np.array([GCVscore(p) for p in grid])
            best = np.argmin(scores)
            if best == 0:
                return y.copy()
            if best == len(grid) - 1:
                s = 10**grid[best]
            else:
                result = minimize_scalar(GCVscore, bounds=(grid[best-1], grid[best+1]), method='bounded')
                if not result.success:
                    raise RuntimeError("GCV optimization failed.")
                s = 10**result.x
            Gamma = 1 / (1 + s * Lambda_sq)
            new_z = idct(Gamma * DCTy, norm='ortho')
            tol = np.linalg.norm(new_z - z) / max(np.linalg.norm(new_z), 1.e-15)
            z = new_z
            if tol <= 1.e-5:
                break
        else:
            raise RuntimeError("Smoothing iteration failed to converge.")

        residual = data - z
        MAD = np.median(np.abs(residual - np.median(residual)))
        if MAD <= 1.e-12:
            break
        tmp = np.sqrt(1 + 16 * s)
        h = np.sqrt(1 + tmp) / np.sqrt(2) / tmp
        u = np.abs(residual) / (1.4826 * MAD * np.sqrt(max(1 - h, 1.e-15))) / 4.685
        W = np.where(u < 1, (1 - u**2)**2, 0)
        if np.sum(W) == 0:
            raise RuntimeError("All smoothing weights vanished.")

    return center + scale * z