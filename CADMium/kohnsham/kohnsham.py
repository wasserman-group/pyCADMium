"""
kohnsham.py
"""

from dataclasses import dataclass
from multiprocessing import Process, Manager, current_process
from pydantic import validator, BaseModel
import numpy as np

from .scf import scf
from ..common.coulomb import coulomb
from ..libxc.libxc import Libxc
from ..hartree.hartree import Hartree
from ..pssolver.pssolver import Pssolver

@dataclass
class V:
    pass

@dataclass
class E:
    E  : float = 0.0 
    Ec : float = 0.0
    Ex : float = 0.0
    Eks: np.ndarray = np.empty((1,1))
    Vks: np.ndarray = np.empty((1,1))

class KohnShamOptions(BaseModel):
    interaction_type : str = 'dft'
    sym : bool = False
    fractional : bool = False
    xfunc_id : int = 1
    cfunc_id : int = 12
    xc_family : str = 'lda'

    @validator('interaction_type')
    def interaction_type_values(cls, v):
        values = ['dft', 'ni']
        if v not in values:
            raise ValueError(f"'interaction_type' must be one of the options: {values}")
        return v

class Kohnsham():
    """
    Handles a standard Kohn-Sham calculation
    """

    def __init__(self, grid, Za, Zb, pol, Nmo, N, optKS={}):

        #Validate options
        optKS =  {k.lower(): v for k, v in optKS.items()}
        for i in optKS.keys():
            if i not in KohnShamOptions().dict().keys():
                raise ValueError(f"{i} is not a valid option for KohnSham")
        optKS = KohnShamOptions(**optKS)

        #Options
        self.optKS = optKS

        #Coordinates
        self.grid = grid

        self.Nmo = np.array(Nmo)
        self.N = np.array(N)

        assert self.Nmo.shape == self.N.shape, "Nmo should be same size as N"

        #Polarization for fragments
        self.pol = pol
        
        #Structures to store component potentials and energies
        self.V = V()
        self.E = E()

        #Handle for the solver object
        self.solver = Pssolver(self.grid, self.Nmo, self.N,
                             {"fractional" : optKS.fractional, "sym" : optKS.sym})

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
        self.Q = np.zeros((self.grid.Nelem))
        self.Alpha = None
        self.Beta = None

        #Calculates coulomb potential corresponding to Za and Zb
        self.calc_nuclear_potential()

        #Libxc/Hartree handles for fragment calculations
        if optKS.interaction_type == 'dft':
            self.exchange = Libxc(self.grid, optKS.xc_family, optKS.xfunc_id)
            self.correlation = Libxc(self.grid, optKS.xc_family, optKS.cfunc_id)
            self.hartree = Hartree(grid)
        elif optKS.interaction_type == 'hfs':
            self.exchange = Libxc(self.grid, optKS.xc_family, optKS.xfunc_id)
            self.correlation = 0.0
            self.hartree = Hartree(grid)
        else:
            self.exchange = 0.0
            self.correlation = 0.0
            self.hartree = 0.0

        #Loop through array and setup solver objects
        for i in range(self.Nmo.shape[0]):
            for j in range(self.Nmo.shape[1]):
                self.solver[i,j].hamiltonian()

                if self.optKS.interaction_type == "ni":
                    self.solver[i,j].e0 = -1.5 * max(self.Za, self.Zb)**2 / (self.solver[i,j].m + 1)**2
                else:
                    self.solver[i,j].e0 = - max(self.Za, self.Zb)**2 / (self.solver[i,j].m + 1)**2 


    # def scf(self, optKS):
    #     scf(self, optKS)
    def scf(self, ks_scf_options={}):
        scf(self, ks_scf_options)

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

        if self.optKS.sym is True:
            self.vnuc = 0.5 * (self.vnuc + self.grid.mirror(self.vnuc))
        
    def set_effective_potential(self):
        """
        Sets new effective potential
        """

        self.veff = self.vnuc + self.vext + self.vhxc
        for j in range(self.Nmo.shape[1]):
            for i in range(self.Nmo.shape[0]):
                self.solver[i,j].setveff(self.veff[:,j])

    def set_veff_external(self, veff):
        """
        Sets user defined veff to solvers
        """
        for j in range(self.Nmo.shape[1]):
            for i in range(self.Nmo.shape[0]):
                self.solver[i,j].setveff(self.veff[:,j])

    def calc_density(self, ITERATIVE=False, dif=0.0):
        #Removed setting Iterative False if only one argument is given

        starttol = 0.01

        #Calculate new densities      
        nout = np.zeros((self.grid.Nelem, self.pol))  

        processes   = []
        manager     = Manager()
        eig_results = manager.dict() 

        for j in range(self.Nmo.shape[1]):
            for i in range(self.Nmo.shape[0]):
                if ITERATIVE is True and dif < starttol:
                    process = Process(target=self.solver[i,j].iter_orbitals,
                                       args=((i,j), eig_results))  
                    processes.append(process)
                    process.start()
                else:
                    process = Process(target=self.solver[i,j].calc_orbitals,
                                      args=((i,j), eig_results)) 
                    processes.append(process)
                    process.start()

        #Wait for all processes to end
        for process in processes:
            process.join()

        # #Give each solvers their results
        for j in range(self.Nmo.shape[1]):
            for i in range(self.Nmo.shape[0]):
                self.solver[i,j].eig = eig_results[(i,j)][0]
                self.solver[i,j].phi = eig_results[(i,j)][1]

        for j in range(self.Nmo.shape[1]):
            for i in range(self.Nmo.shape[0]):     
                #Calculate Orbital's densities
                self.solver[i,j].calc_density()
                #Add up orbital's densities
                nout[:,j] += np.squeeze(self.solver[i,j].n)

        return nout

    def energy(self):
        #Collect energy
        self.E.Eks = np.zeros((1, self.pol))
        #Collect total potential energy 
        self.E.Vks =  np.zeros((1, self.pol))
        #Collect total kinetic energy
        self.E.Ts  = 0.0
        #Collect all eigenvalues
        self.E.evals = np.empty((0))
    
        #Get KohSham energies from solver object
        for i in range(self.Nmo.shape[0]):
            for j in range(self.Nmo.shape[1]):
                self.solver[i,j].calc_energy()

                self.E.Eks[0,j] += self.solver[i,j].eks
                self.E.Ts       += self.solver[i,j].Ts
                self.E.Vks[0,j] += self.solver[i,j].Vs
                self.E.evals = np.append(self.E.evals, self.solver[i,j].eig)

        if self.optKS.interaction_type == 'ni':
            self.E.Vnuc = self.grid.integrate(np.sum(self.n * self.vnuc, axis=1))
            self.E.Vext = self.grid.integrate(np.sum(self.n * self.vext, axis=1))
            self.E.Enuc = self.grid.integrate(np.sum(self.n * self.vnuc, axis=1))
            self.E.Et = self.E.Ts + self.E.Enuc

        elif self.optKS.interaction_type == 'dft':

            self.E.Enuc = self.grid.integrate(np.sum(self.n * self.vnuc, axis=1))
            self.E.Vext = self.grid.integrate(np.sum(self.n * self.vext, axis=1))
            self.E.Vhxc = self.grid.integrate(np.sum(self.n * self.vhxc, axis=1))

            self.V.vh = self.hartree.v_hartree(self.n)
            self.V.eh = 0.5 * self.V.vh[:,0]
            self.E.Eh = self.grid.integrate(np.sum(self.n, axis=1) * self.V.eh)

            self.E.Ex, self.V.ex, self.V.vx = self.exchange.get_xc(self.n, return_epsilon=True)
            self.E.Ec, self.V.ec, self.V.vc = self.correlation.get_xc(self.n, return_epsilon=True)

            self.E.Et = self.E.Ts + self.E.Enuc + self.E.Eh + self.E.Ex + self.E.Ec

        elif self.optKS.interaction_type == "hfs":

            self.E.Enuc = self.grid.integrate(np.sum(self.n * self.vnuc, axis=1))
            self.E.Vext = self.grid.integrate(np.sum(self.n * self.vext, axis=1))
            self.E.Vhxc = self.grid.integrate(np.sum(self.n * self.vhxc, axis=1))

            self.V.vh = self.hartree.v_hartree(self.n)
            self.V.eh = 0.5 * self.V.vh[:,0]
            self.E.Eh = self.grid.integrate(np.sum(self.n, axis=1) * self.V.eh)

            self.E.Ex, self.V.ex, self.V.vx = self.exchange.get_xc(self.n, return_epsilon=True)
            self.E.Ec, self.V.ec, self.V.vc = 0.0, 0.0, 0.0

            self.E.Et = self.E.Ts + self.E.Enuc + self.E.Eh + self.E.Ex + self.E.Ec


        #Nuclear-Nuclear repulsion
        self.E.Vnn = self.Za * self.Zb / (2 * self.grid.a)
        #Total electronic + nuclear energy
        self.E.E = self.E.Et + self.E.Vnn


        # print("Kinetic Energy", self.E.Ts)
        # print("Nuclear Energy", self.E.Enuc)
        # print("Hartree Energy", self.E.Eh)
        # print("Hartree exchange correlation", self.E.Vhxc)
        # print("Correlation Energy", self.E.Ec)
        # print("Exchange Energy", self.E.Ex)
        # print("External", self.E.Vext)
        # print("Total", self.E.E)
        
    def calc_hxc_potential(self):
        "Calculate new potential for each fragment"

        #Fragment effective potentials

        if self.optKS.interaction_type == "ni":
            self.vhxc = np.zeros_like(self.vext)

        elif self.optKS.interaction_type == "dft":
            #Get effective potential for each element in ensemble
            
            self.E.Ex, self.V.vx = self.exchange.get_xc(self.n)
            self.E.Ec, self.V.vc = self.correlation.get_xc(self.n)
            self.V.vh = self.hartree.v_hartree(self.n)

            self.vhxc = self.V.vx + self.V.vc + self.V.vh

        if self.optKS.sym is True:
            self.vhxc = 0.5 * (self.vhxc + self.grid.mirror(self.vhxc))

    def calc_chempot(self):

        homos = []
        for i in range(self.Nmo.shape[0]):
            for j in range(self.Nmo.shape[1]):
                self.solver[i,j].get_homo()
                homos.append(self.solver[i,j].homo)

        self.u = max(homos)
