"""

finite_difference_2d.py
Build Finite Difference Operator Matrices

"""

from ..common.finite_difference_coefficients import finite_difference_coefficients
import numpy as np

def finite_difference_2d(self):

    #Location
    id = self.i2
    ds = finite_difference_coefficients(id, 2)
    self.d2 = ds
    Nd = len(ds)

    #Set up boundary conditions
    bc = np.zeros((self.bcN, self.bcN))
    for i in range(self.bcN):
        bc += ds[Nd-i-1] * np.diag(np.ones((self.bcN - i)), k=-i)

    Da2 = np.zeros((self.Na, self.Na))
    Dr2 = np.zeros((self.Nr, self.Nr))

    for i in range(Nd):
        Da2 += ds[i] * np.diag(np.ones(self.Na - np.abs(id[i])) , k=id[i])
        Dr2 += ds[i] * np.diag(np.ones(self.Nr - np.abs(id[i])) , k=id[i])

    #Odd Symmetry 
    self.oDa2 = Da2.copy()
    self.oDr2 = Dr2.copy()
    #Add in boundary conditions
    self.oDa2[-1 -self.bcN +1:, -1 - self.bcN +1:] -= np.fliplr(bc)
    self.oDa2[:self.bcN, :self.bcN] -= np.fliplr(bc.T)
    self.oDa2 /= self.ha ** 2
    #Add in boundary conditions
    self.oDr2[:self.bcN, :self.bcN] -= np.fliplr(bc.T)
    self.oDr2 /= self.hr ** 2

    #Even Symmetry
    self.eDa2 = Da2.copy()
    self.eDr2 = Dr2.copy()
    #Add in boundary conditions
    self.eDa2[-1-self.bcN+1:, -1-self.bcN+1:] += np.fliplr(bc)
    self.eDa2[:self.bcN, :self.bcN] += np.fliplr(bc.T)
    self.eDa2 /= self.ha ** 2
    #Add in boundary conditions
    self.eDr2[:self.bcN, :self.bcN] += np.fliplr(bc.T)
    self.eDr2 /= self.hr ** 2
    self.bc2 = bc/self.hr **2

    #Symmetrize
    self.eDa2 = 0.5 * (self.eDa2 + self.eDa2.T)
    self.oDa2 = 0.5 * (self.oDa2 + self.oDa2.T)
    self.eDr2 = 0.5 * (self.eDr2 + self.eDr2.T)
    self.oDr2 = 0.5 * (self.oDr2 + self.oDr2.T)



