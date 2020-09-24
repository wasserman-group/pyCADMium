"""
pssolver.py
"""

import numpy as np
from scipy.sparse import spdiags

from .calc_orbitals import calc_orbitals
from .calc_density import calc_density
from .calc_energy import calc_energy
from .iter_orbitals import iter_orbitals
from .normalize_orbitals import normalize_orbitals
from .get_homo import get_homo


#from .hamiltonian import hamiltonian

eps = np.finfo(float).eps


def Pssolver(grid, Nmo, N,
             FRACTIONAL = False, 
             SYM        = False):

    """
    Generates an array of solvers of shape Nmo[0]xNmo[1]
    """

    Nmo = np.array(Nmo)
    N   = np.array(N)
    solver = np.empty((Nmo.shape[0], Nmo.shape[1]), dtype=object)

    for i in range(Nmo.shape[0]):
        for j in range(Nmo.shape[1]):
            solver[i,j] = i_solver(grid, Nmo[i,j], N[i,j], 
                                i, Nmo.shape[1],
                                FRACTIONAL, SYM) 

    return solver
    

class i_solver():
    """
    Keeps track of eigenvalues/vectors and constructs densities and responses.
    """
    def __init__(self, grid, Nmo, N, 
                             Nlvls, pol,
                FRACTIONAL = False,
                SYM        = False):

        verbose=False

        self.grid = grid
        self.Nmo = Nmo
        self.N = N

        #Polarization of electrons handled by this solver
        self.pol = pol

        #Effective Potential
        self.veff = np.zeros((self.grid.Nelem, self.pol))

        #Base Hamiltonian
        self.H0 = None

        #Molecular Orbitals
        self.phi = None
        #Eigenvalues
        self.eig = None
        
        #Estimate of lowest energy value
        self.e0 = None
        
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

        self.opt = {"tol" : 1e-16, "v0":np.ones((self.grid.Nelem, 1))}

        self.Nlvls = Nlvls
        self.pol = pol

        if self.N == -1:
            if self.polarization == 1:
                self.N = 2 * self.Nmo
            else:
                self.N = self.Nmo

        if verbose is True:
            print("\n Warning: Polarization from PSsolver may not ready. Check 'Fill in default number of electrons'")

        self.m = Nlvls

    def hamiltonian(self):
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
            #Odd symmetry
            oT = -0.5 * self.grid.olap
            self.H0 = oT + self.m ** 2 * W @ f 

    def calc_orbitals(self):
        calc_orbitals(self)

    def iter_orbitals(self):
        iter_orbitals(self)

    def normalize_orbitals(self):
        normalize_orbitals(self)

    def calc_density(self):
        calc_density(self)

    def calc_energy(self):
        calc_energy(self)

    def get_homo(self):
        get_homo(self)

    def setveff(self, veff):
        """
        Distribute effective potential to solver object
        Size of Potential array should match polarization
        Assert solver(0).pol == len(veff, 1)
        """
        self.veff = veff