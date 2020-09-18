"""
calc_energy.py
"""

import numpy as np

def calc_energy(self):

    #Initialize energy
    self.eks = 0.0

    #Figure out number of occupied orbitals : Nocc
    if self.m == 0:
        if self.pol == 1: #Unpolarized electrons
            Nocc = np.floor(self.N / 2)
            nu = self.N / (2 - Nocc)

        else: #Polarized electrons
            Nocc = np.floor(self.N)
            nu = self.N - Nocc

    else: # m>0 orbitals hold twice as many elctrons due to plus minus symmetry
        if self.pol == 1: #Unpolarized electrons
            Nocc = np.floor(self.N / 4)
            nu = self.N / (4 - Nocc)

        else: #Polarized electrons
            Nocc = np.floor(self.N / 2)
            nu = self.N / (2.0 - Nocc)

    #Sum energy
    for i in range(int(Nocc)):
        self.eks += self.eig[i]

    #If we are doing fractional orbitals and are non-integer:
    if self.FRACTIONAL is True and nu != 0:
        self.eks += nu * self.eig[int(Nocc)]

    #Scale energy appropriately 
    if self.m == 0:
        if self.pol == 1:
            self.eks = 2 * self.eks

    else:
        if self.pol == 1:
            self.eks = 4 * self.eks
        else:
            self.eks = 2 * self.eks



    #Calculate potential energy of Kohn Sham 
    self.calc_density()
    self.Vs = self.grid.integrate((self.veff[None].T * self.n)[:,0])
    self.Ts = self.eks - self.Vs

    if np.isnan(self.Vs) is True:
        self.Vs = 0.0