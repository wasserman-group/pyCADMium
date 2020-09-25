"""
linresponse.py
"""

import sys

import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import spdiags

def Ws(vs, self):
    """
    Calculates G and its gradient for a given potential
    """

    if self.optInversion["AB_SYM"] is True:
        vs = 0.5 * (vs + self.grid.mirror(vs))

    #Transfer new potentials to solver objects and calculate new densities
    self.solver[0,0].setveff(vs)
    self.solver[0,0].calc_orbitals()
    self.solver[0,0].calc_density()
    self.solver[0,0].calc_energy()
    self.solver[0,0].calc_response()

    #Calculate new density      
    n = np.zeros((self.grid.Nelem, self.pol))  
    for i in range(self.Nmo.shape[0]):
        for j in range(self.Nmo.shape[1]):
            n[:,j] += np.squeeze(self.solver[i,j].n)

    if self.optInversion["AB_SYM"] is True:
        n = 0.5 * (n + self.grid.mirror(n))

    #Calculate error function
    grad = np.hstack( ( B @ (n-n0), self.vs ) )

    #Calculate jacobian of error function
    Jac     = np.hstack( ( B @ self.solver[0,0].chi, np.zeros((1, Nelem)) ) )
    Jac[-1] = 1

    if optInversion["AB_SYM"] is True:
        Jac[-1, :] = 0.5 * (Jac[-1,:] + self.grid.mirror(Jac[-1,:]))

    return grad

def linresponse(self, n0, vs0=None):
    """
    wuyang like inversion of the density with response
    """

    n0 = n0[:, None]
    pol   = n0.shape[1] if len(n0.shape) > 1 else 1
    Nelem = n0.shape[0]
    w = self.grid.w if pol == 1 else np.hstack((self.grid.w, self.grid.w))
    
    if self.solver[0,0].veff is None:
        if vs0 == None:
            #If no initial guess is provided
            #Use the von weizacker inversion
            vs0 = (0.5 * self.grid.elap @ (n0**0.5)) / (n0**0.5 * w) 
            vs0 -= vs[-1]
    else:
        vs0 = self.solver[0,0].veff

    self.vs = np.zeros_like(vs0)
    self.us = np.zeros((1, pol))

    if self.optInversion["ENS_SPIN_SYM"] is True:
        #Invert density n0 to find vs
        i = 1
        #Preallocation
        #B is the inverse of n0 in main diagonal
        B = spdiags( 1./n0[:,i],0, Nele, Nelem)

        self.solver[0,0].hamiltonian()
        self.solver[0,0].e0 = -20
        self.vs[:, i] = least_squares(Ws, vs0[:, i])

        self.us = self.solver[0,0].get_homo()
        self.vs[:,1] = self.vs[:, 0]
        self.us[0,1] = self.us[0,0]
        self.solver[0,0].setveff(self.vs[:, 1])

    else:
        for i in range(pol):
            #Invert density n0 to find vs

            #Preallocation
            B = spdiags(1./n0[:, i], 0, Nelem, Nelem)

            self.solver[0,0].hamiltonian()
            self.solver[0,0].e0 = -20
            self.vs[:, i] = least_squares(Ws, vs0[:,i])
            self.us = self.solver[0,0].get_homo()





        

        

