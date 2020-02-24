"""

NC.py
Calculate Newton_Cotes weights for 1d integration

"""

import numpy as np

def NC(N, Nm):

    xi = np.array(range(N))
    n_p = np.array(range(N))
    x, npow = np.meshgrid(xi, n_p)
    A = x**npow
    b = 1./(n_p+1) * (N-1)**(n_p+1)
    w = np.linalg.lstsq(A,b, rcond=None)
    w = w[0]

    xi = np.array(range(N), dtype=float)
    xi[1:] = xi[1:] - 0.5
    n_p = np.array(range(N))
    x, npow = np.meshgrid(xi, n_p)
    A = x**npow
    b =  1./(n_p+1) * (.5)**(n_p+1)
    w1 = np.linalg.lstsq(A,b, rcond=None)
    w1 = w1[0][1:]


    NP = (N-1) * Nm + 1
    wi = np.zeros(NP)
    wi[:N-1] = w1[:]
    wi[-(N-1):] = w1[::-1]

    for i in range(Nm):
        wi[np.array(range(N))+(N-1)*(i)] = wi[np.array(range(N))+(N-1)*(i)] + w

    return wi