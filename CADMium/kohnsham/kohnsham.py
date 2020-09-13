"""
kohnsham.py
"""

import numpy as np

from .scf import scf
from ..common.coulomb import coulomb
from ..libxc.libxc import Libxc
from ..hartree.hartree import Hartree
from ..pssolver.pssolver import Pssolver

class V:
    pass

class E:
    pass

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

        self.Nmo = np.array(Nmo)
        self.N = np.array(N)

        assert self.Nmo.shape == self.N.shape, "Nmo should be same size as N"

        #Polarization for fragments
        self.pol = pol
        
        #Structures to store component potentials and energies
        self.V = V
        self.E = E

        #Handle for the solver object
        self.solver = np.empty((self.Nmo.shape[0], self.Nmo.shape[1]), dtype=object)
        
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
        for i in range(self.Nmo.shape[0]):
            for j in range(self.Nmo.shape[1]):
                i_solver = Pssolver(self.grid, self.Nmo[i,j], self.N[i,j], 
                                    i, self.Nmo.shape[1],
                                    self.optKS["FRACTIONAL"], self.optKS["SYM"])
                i_solver.hamiltionian()
                self.solver[i,j] = i_solver

        for i in range(self.Nmo.shape[0]):
            for j in range(self.Nmo.shape[1]):
                if self.optKS["interaction_type"] == "ni":
                    self.solver[i,j].e0 = -1.5 * max(self.Za, self.Zb)**2 / (self.solver[i,j].m + 1)**2
                else:
                    self.solver[i,j].e0 = - max(self.Za, self.Zb)**2 / (self.solver[i,j].m + 1)**2 


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
        for i in range(self.Nmo.shape[0]):
            for j in range(self.Nmo.shape[1]):
                self.solver[i,j].setveff(self.veff)

    def calc_density(self, ITERATIVE=False, dif=0.0):
        #Removed setting Iterative False if only one argument is given

        starttol = 0.01

        #Calculate new densities      
        nout = np.zeros((self.grid.Nelem, self.pol))  
        for i in range(self.Nmo.shape[0]):
            for j in range(self.Nmo.shape[1]):
                if ITERATIVE is True and dif < starttol:
                    self.solver[i,j].iter_orbitals()
                else:
                    self.solver[i,j].calc_orbitals()

                #Calculate Orbital's densities
                self.solver[i,j].calc_density()

                #Add up orbital's densities
                nout[:,j] += np.squeeze(self.solver[i,j].n)

        return nout

    def calc_energies(self):
        
        #Get KohSham energies from solver object
        for i in range(self.Nmo.shape[0]):
            for j in range(self.Nmo.shape[1]):
                self.solver[i,j].calc_energy()
                self.E.Eks = self.solver[i,j].get_eks()
                self.E.Ts = self.sovler[i,j].get_Ts()
                self.E.Vks = self.solver[i,j].get_Vs()

        self.E.evals =
        

    def calc_chempot(self):

        homos = []
        for i in range(self.Nmo.shape[0]):
            for j in range(self.Nmo.shape[1]):
                self.solver.[i,j].get_homo()
                homos.append(self.solver[i,j].homo)

        self.u = max(homos)

    def calc_hxc_potential(self):
        "Calculate new potential for each fragment"

        #Fragment effective potentials

        if self.optKS["interaction_type"] == "ni":
            self.vhxc = np.zeros_like(self.vext)

        elif self.optKS["interaction_type"] == "dft":
            #Get effective potential for each element in ensemble
            _, self.V.vx = self.exchange.get_xc(self.n)
            _, self.V.vc = self.correlation.get_xc(self.n)
            self.V.vh = self.hartree.v_hartree(self.n)

            self.vxc = self.V.vx + self.V.vc + self.V.vh

        if self.optKS["SYM"] is True:
            self.vhxc = 0.5 * (self.vhxc + self.grid.mirror(self.vxc))





        



        
        


            

        
         