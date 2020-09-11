"""
kohnsham.py
"""

import numpy as np

from .scf import scf
from ..common.coulomb import coulomb
from ..libxc.libxc import Libxc
from ..hartree.hartree import Hartree
from ..pssolver.pssolver import Pssolver

class Kohnsham():
    """
    Handles a standard kohn sham calculation
    """

    def __init__(self, grid, Za, Zb, pol, Nmo, N, optKS):

        optKS["interaction_type"] = optKS["interaction_type"] if "interaction_type" in optKS.keys() else "dft"
        optKS["SYM"] = optKS["SYM"] if "SYM" in optKS.keys() else False
        optKS["FRACTIONAL"] = optKS["FRACTIONAL"] if "FRACTIONAL" in optKS.keys() else False
        optKS["x_func_id"] = optKS["x_func_id"] if "x_func_id" in optKS.keys() else 1
        optKS["c_func_id"] = optKS["c_func_id"] if "c_func_id" in optKS.keys() else 12
        optKS["xc_family"] = optKS["xc_family"] if "xc_family" in optKS.keys() else "lda"

        #Options
        self.optKS = optKS

        #Coordinates
        self.grid = grid

        #Calculation types / Fragment specification
        #self.interaction_type = optKS["interaction_type"]

        #DFT options for fragment calculations
        #self.xc_family = optKS["xc_family"]
        #self.x_func_id = optKS["x_func_id"]
        #self.c_func_id = optKS["c_func_id"]

        self.Nmo = Nmo
        self.N = N

        assert self.Nmo.shape == self.N.shape, "Nmo should be same size as N"

        #Polarization for fragments
        self.pol = pol
        
        #Structures to store component potentials and energies
        self.V = None
        self.E = None

        #Handle for the solver object
        self.solver = []
        
        #Nuclear charges
        self.Za = Za
        self.Zb = Zb

        #Potentials
        self.vnuc = None
        self.vext = None
        self.vhxc = None
        self.veff = None

        #Density
        self.n = None
        #Chemical potential
        self.u = None 

        #Scaling factors for use with ensembles
        #Fragment energies/densities are scaled by this amount
        self.scale = None
        #Qfunction
        self.Q = None

        #Flags
        #Use AB symmetry for homonuclear diatomics
        # If True, potentails and densities will eb symmetrized
        self.optKS["SYM"] = False
        #Allow fractional occupation of the Homo
        self.optKS["FRACTIONAL"] = False

        self.Alpha = None
        self.Beta = None

        #Calculates coulomb potential corresponding to Za and Zb
        self.calc_nuclear_potential()

        #Libxc/Hartree handles for fragment calculations
        if self.optKS["interaction_type"] == 'dft':
            self.exchange = Libxc(self.grid, self.optKS["xc_family"], self.optKS["x_func_id"])
            self.correlation = Libxc(self.grid, self.optKS["xc_family"], self.optKS["x_func_id"])
            self.hartree = Hartree(grid)
        else:
            self.exchange = 0.0
            self.correlation = 0.0
            self.hartree = 0.0

        #Loop through array and setup solver objects
        for i in range(self.Nmo.shape[1]):
            i_solver = Pssolver(self.grid, self.Nmo, self.N, self.optKS["FRACTIONAL"], self.optKS["SYM"])
            i_solver.hamiltionian()
            self.solver.append(i_solver)
        
        for i_solver in self.solver:
            if self.optKS["interaction_type"] == "ni":
                i_solver.e0 = -1.5 * max(self.Za, self.Zb)**2 / (i_solver.m + 1)**2
            else:
                i_solver.e0 = - max(self.Za, self.Zb)**2 / (i_solver.m + 1)**2 

    def scf(self, optKS):
        scf(self, optKS)

    def calc_nuclear_potential(self):
        """
        Calculate nuclear potential
        """

        v  = coulomb(self.grid, self.Za, self.Zb)

        self.vnuc = np.zeros((v.shape[0], self.pol))

        if self.pol == 1:
            self.vnuc[:, 0] = v

        elif self.pol == 2:
            self.vnuc[:, 0] = v
            self.vnuc[:, 1] = v

        if self.optKS["SYM"] is True:
            self.vnuc = 0.5 * (self.vnuc + self.grid.mirror(self.vnuc))
        
    def set_effective_potential(self):
        """
        Sets new effective potential
        """
        self.veff = self.vnuc + self.vext + self.vhxc
        for i in self.solver:
            i.setveff(self.veff)
         