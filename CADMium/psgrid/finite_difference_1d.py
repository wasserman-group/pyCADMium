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
    bc = np.zeros((self.bcN, self.bcN))
    for i in range(self.bcN):
        bc += ds[Nd-i-1] * np.diag(np.ones((self.bcN - i)), k=-i)

    Da1 = np.zeros((self.Na, self.Na))
    Dr1 = np.zeros((self.Nr, self.Nr))

    for i in range(Nd):
        Da1 += ds[i] * np.diag(np.ones(self.Na - np.abs(id[i])) , k=id[i])
        Dr1 += ds[i] * np.diag(np.ones(self.Nr - np.abs(id[i])) , k=id[i])

    #Antisymmetrize
    Da1 = 0.5 * (Da1 - Da1.T)
    Dr1 = 0.5 * (Dr1 - Dr1.T)

    #Odd symmetry
    self.oDa1 = Da1.copy()
    self.oDr1 = Dr1.copy()
    #Add in boundary conditions for angular derivatives
    self.oDa1[-1 - self.bcN + 1:, -1 - self.bcN +1 :] -= np.fliplr(bc) 
    self.oDa1[:self.bcN, :self.bcN] += np.fliplr(bc.T)
    self.oDa1 /= self.ha
    #Add in boundary conditions for radial derivatives
    self.oDr1[:self.bcN, :self.bcN] += np.fliplr(bc.T)
    self.oDr1 /= self.hr

    #Even symmetry 
    self.eDa1 = Da1.copy()
    self.eDr1 = Dr1.copy()
    #Add in boundary conditions for angular derivatives
    self.eDa1[-1 - self.bcN + 1:, -1 - self.bcN +1:] += np.fliplr(bc)
    self.eDa1[: self.bcN, :self.bcN] -= np.fliplr(bc.T)
    self.eDa1 /= self.ha 
    #Add in boundary conditions for radial derivatives
    self.eDr1[:self.bcN, :self.bcN] -= np.fliplr(bc.T)
    self.eDr1 /= self.hr
    
    self.bc1 = bc / self.hr



