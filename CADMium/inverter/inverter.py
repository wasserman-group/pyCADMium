"""
inverter.py
"""

import numpy as np
from pydantic import validator, BaseModel

from .linresponse import linresponse
from .orbitalinvert import orbitalinvert
from .simple import simple

from .get_ts_WFI import get_ts_WFI
from .get_Ts import get_Ts
# from .linresponse import Ws
# from .linresponse import Jacobian

class InverterOptions(BaseModel):
    ab_sym : bool = False
    ens_spin_sym : bool = False
    use_iterative : bool = False
    disp : bool = False
    avoid_loop : bool = False
    invert_type : str = 'wuyang'
    k_val : float = -1e-12
    tol_lin_solver : float = 1e-2
    tol_invert : float = 1e-15
    max_iter_lin_solver : int = 2000
    max_iter_invert : int = 20
    res_factor : float = 1e0

    @validator('invert_type')
    def invert_type_values(cls, v):
        values = ['wuyang', 'simple', 'orbitalinvert', 'qinvert', 'eigensolveinvert', 'test']
        if v not in values:
            raise ValueError(f"'invert_type' must be one of the options: {values}")
        return v

class Inverter():
    """
    Inverter handles inversion algorithms

    Parameters 
    ----------

    Attributes
    ----------
    vs : numpy array
        Kohn Sham Potential
    us : float
        Chemical Potential
    ts_WFI :
        Kinetic energy density using laplacian
    ts_WFII :
        Kinetic energy density using gradient
    Ts : float
        Kinetic Energy


    Defaults
    --------
    Convergence and Algorithm Parameters

    Use_Iterative : logical
        Set to True to use iterative sovler with orbital invert or eigesolveinvert.
        Otherwise direct solver will be used which may take a long time for large
        grids.

    Kval : floag
        Parameter for iterative preconditioner

    TolLinsolve : float
        Tolerance for iterative solution to linear solve

    MaxIterLinsolve: int
        Maximum number of iterations for solution to linear solve

    TolInvert: float
        Requested tolerance for inversion

    MaxIterInvert: int
        Maximum iterations for inversions

    ResFactor: 1e0
        Determines choice of update in orbitalinvert
        
    """

    def __init__(self, grid, solver, optInv={}):

        optInv =  {k.lower(): v for k, v in optInv.items()}
        for i in optInv.keys():
            if i not in InverterOptions().dict().keys():
                raise ValueError(f"{i} is not a valid option for Inverter")
        optInv = InverterOptions(**optInv)
        self.optInv = optInv
    
        self.grid = grid
        self.Nelem = grid.Nelem
        self.solver = solver
        self.Nmo = solver[0,0].Nmo

        self.B  = None
        self.n0 = None

        self.vs = None
        self.us = None
        self.ts_WFI = None
        self.ts_WFII = None

    def invert(self, n0, vs0, phi0=[], e0=[], ispin=[], Qi=[]):
        """
        Do the inverstion
        """

        if self.optInv.invert_type == "wuyang":
            flag, output = self.linresponse(n0, vs0, ispin)

        elif self.optInv.invert_type == "simple":
            flag, output = self.simple(n0, vs0, ispin)

        elif self.optInv.invert_type == "orbitalinvert":
            flag, output = self.orbitalinvert(n0, vs0, phi0, e0, ispin)

        elif self.optInv.invert_type == "qinvert":
            flag, output = self.qinvert(n0, vs0, phi0, e0, ispin, Qi)

        elif self.optInv.invert_type == "eigensolveinvert":
            flag, output = self.eigensolveinvert(n0, vs0, ispin)

        elif self.optInv.invert_type == "test":
            flag, output = self.test(n0, vs0, phi0, e0, ispin)
        else:
            raise ValueError(f"{self.optInv.invert_type} is not an available inversion method")

        return flag, output

    def get_vt(self):
        """
        Gets kinetic potential
        """
        vt = np.ones((self.grid.Nelem, 1)) * self.us - self.vs
        return vt

    def get_Ts(self):
        Ts = get_Ts(self)
        return Ts

    def get_ts_WFI(self):
        ts = get_ts_WFI(self)
        return ts

    def simple(self, n0, vs0, ispin):
        flag, output = simple(self, n0, vs0, ispin)
        return flag, output

    def orbitalinvert(self, n0, vs0, phi0, e0, ispin):
        flag, output = orbitalinvert(self, n0, vs0, phi0, e0, ispin)
        return flag, output

    #Inversion methods
    def linresponse(self, n0, vs0, ispin):
        flag, output = linresponse(self, n0, vs0, )
        return flag, output

    def Ws(self, vs, spin):
        """
        Calculates G for a given potential
        """
        if self.optInv.ab_sym is True:
            vs = 0.5 * (vs + self.grid.mirror(vs))

        # fig = plt.figure()
        # x,y = self.grid.axis_plot(vs)
        # plt.plot(x,y)
        # # plt.xlim(-7,7)
        # # plt.ylim(-5,5)
        # plt.show()

        #Transfer new potentials to solver objects and calculate new densities
        self.solver[0,spin].setveff(vs)
        self.solver[0,spin].calc_orbitals()
        self.solver[0,spin].calc_density()
        self.solver[0,spin].calc_energy()
        self.solver[0,spin].calc_response()

        #Calculate new density     
        n = np.zeros((self.grid.Nelem, 1))  
        for i in range(self.solver.shape[0]):
            n[:,0] += np.squeeze(self.solver[i,spin].n)

        if self.optInv.ab_sym is True:
            n = 0.5 * (n + self.grid.mirror(n))

        if np.isnan(self.vs).any():
            print("vs is nan")
        if np.isinf(self.vs).any():
            print("vs is inf")

        #Calculate error function
        grad = np.vstack( ( self.B @ (n-self.n0[:, spin][:,None]), self.vs[self.Nelem-1, spin] ) )
        grad = np.squeeze(grad)

        if np.isnan(grad).any() is True:
            print("Grad is nan")
        if np.isinf(grad).any() is True:
            print("Grad is inf")

        # print("max grad", np.linalg.norm(grad))

        return grad

    def Jacobian(self, vs, spin):
        """
        Calculates Jacobian for a given vs
        """


        #Calculate jacobian of error function.
        Jac = np.vstack( ( self.B @ self.solver[0,spin].chi, np.zeros((1, self.Nelem)) ) )
        Jac = np.asarray(Jac)
        Jac[-1, -1] = 1

        if self.optInv.ab_sym is True:
            Jac[-1,:] = 0.5 * (Jac[-1,:] + self.grid.mirror(Jac[-1,:]))

        if np.isnan(Jac).any():
            print("Jac is nan")
        if np.isinf(Jac).any():
            print("Jac is inf")

        return Jac
