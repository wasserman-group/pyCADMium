"""
sigma.py
"""

import numpy as np 

def sigma(self, n):

    pol = n.shape[1] 

    if pol == 1:
        sigma = (self.grada @ n) ** 2 + (self.gradr @ n) ** 2

    elif pol == 2:
        npoints = n.shape[0]
        sigma = np.zeros((npoints, 3))
        sigma[:,2] = (self.grada @ n[:, 1]) ** 2 +  (self.gradr @ n[:, 1]) ** 2
        sigma[:,0] = (self.grada @ n[:, 0]) ** 2 +  (self.gradr @ n[:, 0]) ** 2
        sigma[:,1] = self.grada @ n[:, 0] * self.grada @ n[:,1] + \ 
                     self.gradr @ n[:, 0] * self.gradr @ n[:,1]

    else:
        raise ValueError("Shape of density in second axis should not be greater than 2")

    return sigma