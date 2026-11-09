# Coordinate conversions used by the course scripts
import sys
from pathlib import Path
import numpy as np

# Course convention is [a, e, i, w, OM, theta], input angles are radians
sys.path.append(str(Path(__file__).resolve().parents[1] / 'curtis_scripts' / 'python'))
from sv_from_coe import sv_from_coe
from coe_from_sv import coe_from_sv

def orb2xyz(mu, oe):
    a, e, i, w, OM, theta = oe
    p = a*(1 - e**2)
    if p <= 0:
        raise ValueError("These course elements require a positive semilatus rectum.")
    h = np.sqrt(mu*p)
    r_vec, v_vec = sv_from_coe([h, e, OM, i, w, theta], mu)
    return np.concatenate((r_vec, v_vec))

def xyz2orb(mu, r_vec, v_vec):
    h, e, OM, i, w, theta, a = coe_from_sv(r_vec, v_vec, mu)
    # CSV output uses degrees, as indicated by the course-script column labels
    return [a, e, np.rad2deg(i), np.rad2deg(w), np.rad2deg(OM), np.rad2deg(theta)]