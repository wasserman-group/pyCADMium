"""
linresponse.py
"""

import numpy as np
from scipy.optimize import least_squares
#from scipy.optimize import minimize
from scipy.sparse import spdiags


def linresponse(self, n0, vs0=None):
    """
    wuyang like inversion of the density with response
    """

    n0 = n0[:, None] if len(n0.shape) == 1 else n0
    pol   = n0.shape[1] if len(n0.shape) > 1 else 1
    self.pol = pol
    self.n0  = n0
    Nelem = n0.shape[0]
    w = self.grid.w if pol == 1 else np.hstack((self.grid.w, self.grid.w))
    
    if self.solver[0,0].veff is None:
        if vs0 is None:
            #If no initial guess is provided
            #Use the von weizacker inversion
            vs0 = (0.5 * self.grid.elap @ (n0**0.5)) / (n0**0.5 * w) 
            vs0 -= vs[-1]
    else:
        vs0 = self.solver[0,0].veff[:,None] if len(self.solver[0,0].veff.shape) == 1 else self.solver[0,0].veff

    self.vs = np.zeros_like(vs0)
    self.us = np.zeros((1, pol))
    flag   = np.empty_like(self.solver, dtype=object)
    output = np.empty_like(self.solver, dtype=object)
    
    if self.optInv.ens_spin_sym is True:
        #Invert density n0 to find vs
        #Preallocation
        #B is the inverse of n0 in main diagonal
        B = spdiags( 1./n0[:,0],0, Nelem, Nelem)

        self.B = B

        for j in range(self.solver.shape[0]):
            self.solver[j,0].hamiltonian()
            self.solver[j,0].e0 = -20

        res_lsq = least_squares(fun    = self.Ws,
                                x0     = vs0[:, 0], 
                                jac    = self.Jacobian, 
                                method = "trf", 
                                args   = (0,),
                                xtol   = self.optInv.tol_invert,
                                gtol   = self.optInv.tol_invert)

        #Get solution from least squares object
        self.vs[:,0] = res_lsq.x
        self.us[0] = self.solver[0,0].get_homo()
        flag[0,0] = res_lsq.status
        output[0,0] = res_lsq

        #Copy information
        self.vs = np.hstack((self.vs, self.vs))
        self.us = np.hstack((self.us, self.us))
        self.solver[0,1].setveff(self.vs[:,0])
        flag[0,1] = res_lsq.status
        output[0,1] = res_lsq

    else:
        for i in range(pol):
            #Invert density n0 to find vs
            #Preallocation

            B = spdiags( 1./n0[:,i], 0, Nelem, Nelem)
            self.B = B
            for j in range(self.solver.shape[0]):
                self.solver[j,i].hamiltonian()
                self.solver[j,i].e0 = -20

            res_lsq = least_squares(fun    = self.Ws,
                                    x0     = vs0[:, i], 
                                    jac    = self.Jacobian, 
                                    method = "trf", 
                                    args   = (i,),
                                    xtol   = self.optInv.tol_invert,
                                    gtol   = self.optInv.tol_invert)

            for j in range(self.solver.shape[0]):
                flag[j,i] = res_lsq.status
                output[j,i] = res_lsq

            self.vs[:,i] = res_lsq.x
            self.us[i]   = self.solver[0,i].get_homo()


    return flag, output





        

        

