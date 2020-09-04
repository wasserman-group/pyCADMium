"""
kohnsham.py
"""

import numpy as np

from ..common.coulomb import coulomb
from ..pssolver import pssolver

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
        self.x_func_id = partiiton.x_func_id
        self.c_func_id = partition.c_func_id

        if spin == "alpha":
            self.Nmo = partition.Nmo_a
            self.N = partition.N_a
        else:
            self.Nmo = partition.Nmo_b
            self.N = partititon.N_b

        assert(self.Nm0, self.N)

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
        self.solver = None
        
        #Nuclear charges
        self.Za = partition.Za
        self.Zb = partition.Zb

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
        self.SYM = partition.SYM
        #Allow fractional occupation of the Homo
        self.FRACTIONAL = partition.FRACTIONAL

        self.Alpha = None

        #Calculates coulomb potential corresponding to Za and Zb
        self.partition.calc_nuclear_potential

        if self.interaction_type == 'dft':
            self.exchange = partition.exchange
            self.correlation = partition.correlation
            self.hartree = partition.hartree
        else:
            self.exchange = 0.0
            self.correlation = 0.0
            self.hartree = 0.0

        #Loop through array and setup solver objects
            for i in range(self.Nmo):

        for i in range(self.Nmo):
            i_solver = pssolver(grid, Nmo[i], N[i], self.FRACTIONAL, self.SYM)
            i_solver.hamiltonian()
            self.solver.append(i_solver)

        for i_solver in self.solver:
            if self.interaction_type == "ni":
                i_solver.e0 = -1.5 * np.max(self.Za, self.Zb)**2 / (i_solver.m + 1)**2
            else:
                i_solver.e0 = - np.max(self.Za, self.Zb)**2 / (i_solver.m + 1)**2 


         