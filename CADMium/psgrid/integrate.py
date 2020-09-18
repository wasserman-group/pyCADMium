"""
integrate.py
"""

import numpy as np

def integrate(self, f):
    """
    integrate a function f 
    """

    if len(f.shape) == 1:
        integral = (2 * np.pi * self.hr * self.ha) * np.sum(self.wi * self.w * f)

    elif f.shape[1] == 1:
        integral = (2 * np.pi * self.hr * self.ha) * np.sum(self.wi * self.w * f)

    else:
        integral = np.zeros_like(f)
        for i in range(f.shape[1]):
            integral[:, i] = (2 * np.pi * self.hr * self.ha) * np.sum(self.wi * self.w * f[:, i])

    return integral
