"""
hartree.py
"""

import numpy as np
from scipy.special import lpmv as legendre
from scipy.sparse.linalg import spsolve

class Hartree():
    """
    Handles Calculation of all potentials: Coulomb, Hartree, and Exchange-Correlation
    """

    def __init__(self, grid, 
                            #**kwargs
                            ):
        self.grid = grid

    def v_hartree(self, nh):
        """
        Calculates the hartree potential using density 'nh'
        """

        if len(nh.shape) == 1:
            pol = 1
        else:
            pol = nh.shape[1]
            nh = np.sum(nh, axis=1)

        #Number of multipoles to calculate. 
        #It appears not to affect performance
        n_multipole = 7

        #Fill in mesh and boundary mesh with 'Z' and rho values
        bZ   = self.grid.a * np.cosh(self.grid.bXr) * np.cos(self.grid.bXa)
        #Check the raised to the minus one in next expression
        brho = self.grid.a * (np.cosh(self.grid.bXr)**2 + np.cos(self.grid.bXa)**2-1)**0.5
        Z    = self.grid.a * np.cosh(self.grid.Xr) * np.cos(self.grid.Xa)
        rho  = self.grid.a * (np.cosh(self.grid.Xr)**2 + np.cos(self.grid.Xa)**2-1)**0.5

        #Calculate Multipoles
        #Zero order gets calculated separately
        #Choose order zero of Legendre functions. Hartree potential has m=0 symmetry
        P0 = legendre(0,0, Z/rho)
        bP0 = legendre(0,0, bZ/brho)
        #Integrate to find kth multipole
        Q0 = self.grid.integrate(P0 * nh)
        #Calculate kth multipole contribution to the hartree potential
        #in the boundary region
        bVh = Q0 * np.reciprocal(brho) * bP0

        #Rest of the multipoles:
        for k in range(1,n_multipole+1):
            Pn = legendre(0, k, Z/rho)
            Qn = self.grid.integrate((rho**k) * Pn*nh)
            bPn = legendre(0,k, bZ/brho)
            bVh += Qn * np.reciprocal(brho**(k+1)) * bPn

        #Calculate source term in Poisson's equation
        b  = -4.0 * np.pi * self.grid.w * nh

        #Use lhs of Poisson's equation in the boundary region to find
        #a corresponding source term wich implements boundary conditions
        #See Kobus et al Comp. Phys. Commun. 98(1996) 346-358
        bQ =  (self.grid.blap @ bVh).reshape(self.grid.Na, self.grid.bcN, order='F')
        
        #Add boudnary term into source term. 
        b = b.reshape(self.grid.Na, self.grid.Nr, order='F')
        b[:, -1-self.grid.bcN+1:] = b[:, -1-self.grid.bcN+1:] - bQ
        b = b.reshape(self.grid.Na * self.grid.Nr, 1, order='F')
        
        #Solve discretized Poisson equation using LU decomposed labplacian
        x  = spsolve(self.grid.L_lap, b)
        vh = np.zeros((x.shape[0], pol))
        vh[:, 0] = spsolve(self.grid.U_lap, x)

        if pol == 2:
            vh[:, 1] = vh[:, 0] 

        return vh

    def e_hartree(self, nh):
        """
        Calculate hartree energy per particle
        """

        vh = self.v_hartree(np.sum(nh, axis=1))
        eh = 0.5 * vh

        return eh
