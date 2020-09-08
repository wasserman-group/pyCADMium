"""
pssolver.py
"""

import numpy as np
from scipy.sparse import spdiags

#from .hamiltonian import hamiltonian

eps = np.finfo(float).eps

class Pssolver():
    """
    Keeps track of eigenvalues/vectors and constructs densities and responses.
    """
    def __init__(self, grid, Nmo, N, FRACTIONAL, SYM):

        verbose=False

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
        self.ked_WFI = None #Use laplacian
        self.kedWFII = None #Use gradient
        

        self.FRACTIONAL = FRACTIONAL
        self.SYM = SYM
        self.ITERLINSOLVE = True
        self.TOL_ORBITAL = 2e-15
        self.TOL_IN_SOLVER = 1e-4
        self.tol = eps
        self.v0 = np.ones(self.grid.Nelem)
        self.default_e0 = -20.0

        self.Nlvls = self.Nmo.shape[0]
        self.pol = self.Nmo.shape[1]

        if -1 in self.N:
            if self.polarization == 1:
                self.N = 2 * self.Nmo
            else:
                self.N = self.Nmo

        if verbose is True:
            print("\n Warning: Polarization from PSsolver may not ready. Check 'Fill in default number of electrons'")

        self.m = self.Nlvls - 1
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
        if np.mod(self.m, 2) == 0:
            #Even symmetry
            eT = -0.5 * self.grid.elap
            self.H0 = eT + self.m ** 2 * W @ f
        else:
            oT = -0.5 * self.grid.olap
            self.H0 = oT + self.m ** 2 * W @ f 