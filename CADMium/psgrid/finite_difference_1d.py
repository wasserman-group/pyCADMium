"""

fd1.py
Build Finite Difference Operator Matrices

"""

from ..common.finite_difference_coefficients import finite_difference_coefficients
import numpy as np


def finite_difference_1d(self):

    #Location
    id = self.i1
    ds = finite_difference_coefficients(id, 1)
    self.d1 = ds
    Nd = len(ds)


    #Set up boundary conditions

    bc = np.zeros((self.bcN))

