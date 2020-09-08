"""
kohnsham.py
"""

import numpy as np

from ..common.coulomb import coulomb
from ..pssolver.pssolver import Pssolver

class KohnSham():
    """
    Handles a standard kohn sham calculation
    """

    def __init__(self, partition, spin):

        #Coordinates
        self.grid = partition.grid

        #Calculation types / Fragment specification
        self.interaction_type = partition.interaction_type

        #DFT options for fragment calculations
        self.xc_family = partition.xc_family
        self.x_func_id = partition.x_func_id
        self.c_func_id = partition.c_func_id

        if spin == "alpha":
            self.Nmo = partition.Nmo_a
            self.N = partition.N_a
        elif spin == "beta":
            self.Nmo = partition.Nmo_b
            self.N = partition.N_b

        # if self.Nmo.shape != self.N.shape:
        #     raise ValueError("Shape of Nmo must match N")
        assert(self.Nmo.shape, self.N.shape), "Nmo should be same size as N"

        #Libxc/Hartree handles for fragment calculations
        self.exchange = partition.exchange
        self.correlation = partition.correlation
        self.hartree = None

        #Polarization for fragments
        self.pol = partition.pol
        
        #Structures to store component potentials and energies
        self.V = None
        self.E = None

        #Handle for the solver object
        self.solver = []
        
        #Nuclear charges
        if spin == "alpha":
            self.Za = partition.Za
            self.Zb = 0
        if spin == "beta":
            self.Zb = partition.Zb
            self.Za = 0

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
        self.SYM = False
        #Allow fractional occupation of the Homo
        self.FRACTIONAL = False

        self.Alpha = None

        #Calculates coulomb potential corresponding to Za and Zb
        self.calc_nuclear_potential()

        if self.interaction_type == 'dft':
            self.exchange = partition.exchange
            self.correlation = partition.correlation
            self.hartree = partition.hartree
        else:
            self.exchange = 0.0
            self.correlation = 0.0
            self.hartree = 0.0

        #Loop through array and setup solver objects
        for i in range(self.Nmo.shape[1]):
            i_solver = Pssolver(self.grid, self.Nmo, self.N, self.FRACTIONAL, self.SYM)
            i_solver.hamiltionian()
            self.solver.append(i_solver)
        
        for i_solver in self.solver:
            if self.interaction_type == "ni":
                i_solver.e0 = -1.5 * max(self.Za, self.Zb)**2 / (i_solver.m + 1)**2
            else:
                i_solver.e0 = - max(self.Za, self.Zb)**2 / (i_solver.m + 1)**2 

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

        if self.SYM is True:
            self.vnuc = 0.5 * (self.vnuc + self.grid.mirror(self.vnuc))
        


         