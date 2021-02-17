"""
pssolver.py
"""

import numpy as np
from scipy.sparse import spdiags
from pydantic import validator, BaseModel

from .calc_orbitals import calc_orbitals
from .calc_density import calc_density
from .calc_energy import calc_energy
from .calc_response import calc_response
from .iter_orbitals import iter_orbitals
from .normalize_orbitals import normalize_orbitals
from .get_homo import get_homo

from .calc_ked_WFI import calc_ked_WFI
from .get_ked_WFI import get_ked_WFI

eps = np.finfo(float).eps

class SolverOptions(BaseModel):
    fractional : bool = False
    sym : bool = False
    iter_lin_solver : bool = True
    tol_orbital : float = 1e-15
    tol_lin_solver : float = 1e-4
    disp : bool = True

def Pssolver(grid, Nmo, N, optSolver={}
            #  FRACTIONAL = False, 
            #  SYM        = False
             ):

    """
    Generates an array of solvers of shape Nmo[0]xNmo[1]
    """

    Nmo = np.array(Nmo)
    N   = np.array(N)
    solver = np.empty((Nmo.shape[0], Nmo.shape[1]), dtype=object)

    for i in range(Nmo.shape[0]):
        for j in range(Nmo.shape[1]):
            solver[i,j] = i_solver(grid, Nmo[i,j], N[i,j], 
                                i, Nmo.shape[1], optSolver) 

    return solver
    

class i_solver():
    """
    Keeps track of eigenvalues/vectors and constructs densities and responses.
    """
    def __init__(self, grid, Nmo, N, 
                             Nlvls, pol, optSolver):

        optSolver =  {k.lower(): v for k, v in optSolver.items()}
        for i in optSolver.keys():
            if i not in SolverOptions().dict().keys():
                raise ValueError(f"{i} is not a valid option for KohnSham")
        optSolver = SolverOptions(**optSolver)
        self.optSolver = optSolver

        self.grid = grid
        self.Nmo = Nmo
        self.N = N

        #Polarization of electrons handled by this solver
        self.pol = pol

        #Effective Potential
        #self.veff = np.zeros((self.grid.Nelem, self.pol))
        self.veff = None

        #Base Hamiltonian
        self.H0 = None

        #Molecular Orbitals
        self.phi = None
        #Eigenvalues
        self.eig = None
        
        #Estimate of lowest energy value
        self.e0 = -20
        
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
        

        # self.FRACTIONAL = FRACTIONAL
        # self.SYM = SYM
        # self.ITERLINSOLVE = True
        # self.TOL_ORBITAL = 2e-15
        # self.TOL_IN_SOLVER = 1e-4
        # self.tol = eps
        self.v0 = np.ones(self.grid.Nelem)
        # self.default_e0 = -20.0

        self.opt = {"tol" : 1e-16, "v0":np.ones((self.grid.Nelem, 1))}

        self.Nlvls = Nlvls
        self.pol = pol

        if self.N == -1:
            if self.polarization == 1:
                self.N = 2 * self.Nmo
            else:
                self.N = self.Nmo

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

    def calc_orbitals(self, solver_id=None, return_dict=None):
        calc_orbitals(self, solver_id, return_dict)

    def calc_response(self):
        calc_response(self)

    def iter_orbitals(self, solver_id=None, return_dict=None):
        iter_orbitals(self, solver_id, return_dict)

    def normalize_orbitals(self):
        normalize_orbitals(self)

    def calc_density(self):
        calc_density(self)

    def calc_energy(self):
        calc_energy(self)

    def get_homo(self):
        homo = get_homo(self)
        return homo

    def get_ked_WFI(self):
        ked = get_ked_WFI(self)

    def calc_ked_WFI(self):
        calc_ked_WFI(self)

    # def get_Ts(self):
    #     """
    #     Get total kinetic energy
    #     """
    #     Ts = 0.0
    #     self.calc_energy()
    #     for i in range(self.solver.shape[0]):
    #         for j in range(self.solver.shape[1]):
    #             Ts += np.sum( self.solver[i,j].Ts )

    def setveff(self, veff):
        """
        Distribute effective potential to solver object
        Size of Potential array should match polarization
        Assert solver(0).pol == len(veff, 1)
        """

        self.veff = veff