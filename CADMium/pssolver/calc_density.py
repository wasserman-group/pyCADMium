"""
calc_density.py
"""

import numpy as np

def calc_density(self):
    """
    Calculate density from orbitals
    """
    self.n = np.zeros((self.grid.Nelem,1))

    #Figure out the number of occupied norbitals: Nocc
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

    #Construct Density
    for i in range(int(Nocc)):
        self.n += (self.phi[:,i]**2)[None].T / self.grid.integrate(self.phi[:,i]**2)

    #If we are doing fractional orbitals and are non-integer
    if self.optSolver.fractional and nu != 0:
    # if self.FRACTIONAL and nu != 0:
        self.n += nu * (self.phi[:, int(Nocc)]**2)[None].T / self.grid.integrate(self.phi[:, int(Nocc)]**2)

    #Scale densities appropriately
    if self.m == 0:
        if self.pol == 1:
            self.n = 2 * self.n

    else:
        if self.pol == 1:
            self.n = 4 * self.n
        else:
            self.n = 2 * self.n