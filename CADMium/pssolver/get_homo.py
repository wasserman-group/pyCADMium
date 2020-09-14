"""
get_homo.py
"""

import numpy as np

def get_homo(self):
    "Find  HOMO eigenvalue"

    #Figure out index of homo eigenvalue
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

    if nu != 0:
        Nocc = Nocc+1

    #Record homo eigenvalue
    if Nocc != 0:
        self.homo = self.eig[Nocc]
    else:
        self.homo = -1 / np.finfo(float).eps
    