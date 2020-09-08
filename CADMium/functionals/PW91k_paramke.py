"""
PW91k_paramke
"""

import numpy as np

def PW91k_paramke(s, k):

    A = k[0] @ s * np.asin(k[4] * s)
    B = (k[1] - k[2]@np.exp(-k[3] @ s**2)) * s**2
    C = k[5] @ s**4

    F = (1+A+B) / (1+A+C)

    return F