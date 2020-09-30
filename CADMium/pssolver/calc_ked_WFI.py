"""
calc_ked_WFI
"""

import numpy as np

def calc_ked_WFI(self):
    """
    Calculate kinetic energy density using laplacian of orbitals
    """

    #Initialize kinetic energy density
    self.ked_WFI = np.zeros( (self.grid.Nelem, 1))

    #Figure out the number of occupied orbitals
    if self.m == 0:
        if self.pol == 1:
            Nocc = np.floor(self.N/2)
            nu = self.N / 2 - Nocc
        else:
            Nocc = np.floor(self.N)
            nu = self.N - Nocc

    else:
        #m>0 orbitals hold twice as many electrons due to +-m symmetry
        if self.pol == 1:
            Nocc = np.floor(self.N / 4)
            nu = self.N / 4 - Nocc
        else:
            Nocc = np.floor(self.N/2)
            nu = self.N / 2 - Nocc

    #Construct density
    for i in range(Nocc):
        #Normalized orbital
        phi_norm = self.phi[:,i] / self.grid.integrate( self.phi[:,i]**2)
        self.ked_WFI += (phi_norm * (self.H0 @ phi_norm)) / self.grid.w

    #Scale densities appropriately
    if self.m == 0:
        if self.pol == 1: #Unpolarized electrons
            self.ked_WFI = 2 * self.ked_WFI

    else: # m>0 orbitals hold twice as many electrons due to +-m symmetry
        if self.pol == 1:
            self.ked_WFI = 4 * self.ked_WFI
        else:
            self.ked_WFI = 2 * self.ked_WFI  