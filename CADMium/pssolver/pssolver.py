"""
pssolver.py
"""

import numpy as np
from scipy.sparse import spdiags

#from .hamiltonian import hamiltonian

eps = np.finfo(float).eps

class Pssolver():
    """
    Builds a 
    """

    def __init__(self, grid, Nmo, N, FRACTIONAL, SYM):

        self.grid = grid
        self.Nmo = Nmo
        self.N = N

        #Effective Potential
        self.veff = None

        #Base Hamiltonian
        self.H0 = None

        #Molecular Orbitals
        self.phi = None
        #Eigenvalues
        self.eig = None
        
        #Estimate of lowest energy value
        self.e0 = None

        #Number of MO to calculate
        self.Nmo = None
        #Number of electrons handled by individual solver
        self.N = None
        #Polarization of electrons handled by this solver
        self.pol = None
        
        #KohnSham energy
        self.eks = None
        #Kohn Sham potential 
        self.Vs = None
        #Kohn Sham kinetic energy
        self.Ts = None
        #Electron density
        self.n = None
        #Change in density from given change in potential
        self.chi = None
        #Highest occupied Molecular Orbital
        self.homo = None

        #kinetic energy densities | Not uniquely defined
        self.ked_WFI #Use laplacian
        self.kedWFII #Use gradient
        

        self.FRACTIONAL = FRACTIONAL
        self.SYM = SYM
        self.ITERLINSOLVE = True
        self.TOL_ORBITAL = 2e-15
        self.TOL_IN_SOLVER = 1e-4
        self.tol = eps
        self.v0 = np.ones(self.grid.Nelem)
        self.default_e0 = -20.0

        self.Nlvls = self.Nmo[0]
        self.pol = self.Nmo[1]

        if self.N == -1:
            if self.polarization == 1:
                self.N = 2 * self.Nmo
            else:
                self.N = self.Nmo

        self.m = self.Nvls - 1
        # self.pssolvers = np.array(np.zeros(self.Nlvls, self.pol))

        # self.results = { "Nmo" : np.zeros((self.Nlvls, self.pol)), 
        #                  "N"   : np.zeros((self.Nlvls, self.pol)), 
        #                  "m"   : np.zeros((self.Nlvls, self.pol)), } 
                        
        # for i in range(self.Nlvls):
        #     for j in range(self.pol):
        #         self.results["Nmo"][i,j] = self.Nmo[i,j]
        #         self.results["N"][i,j] = self.N[i,j]
        #         self.results["m"][i,j] = i-1


    def hamiltionian(self):
        """
        Construct basic Hamiltonian H_0
        Includes effective potential due to angular momentum around bond axis
        """

        Nelem = self.grid.Nelem
        
        #Inverse volume element and angular momentum potential
        W = spdiags(data=self.grid.w, diags=0, m=Nelem, n=Nelem)
        f = spdiags(data=self.grid.f, diags=0, m=Nelem, n=Nelem)

        #Choose the correct symmetry for m
        if mod(self.m, 2) == 0
            #Even symmetry
            eT = -0.5 * self.grid
            self.H0 = eT + self.m ** w @ W @ f
        else:
            oT = -0.5 * self.grid.olap
            self.H0 = oT + self.m ** 2 @ W @ f 