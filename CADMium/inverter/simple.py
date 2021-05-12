"""
simple.py
"""

import numpy as np
from scipy.sparse import spdiags


def simple(self, n0, vs0, tilde=None):
    """
    Simple inversion of the density
    Invert Density n0 to vind vs
    """

    pol   = 1 if len(n0.shape) == 1 else 2
    Nelem = n0.shape[0]
    n0 = n0[:None] if len(n0.shape) == 1 else n[:,0][:,None]

    if self.solver[0,0].veff is None:
        vs0 = 0.5 * (0.5 * self.grid.elap * (n0 ** 0.5)) / ( n0**0.5 * self.grid.w)
        vs0 -= vs0[-1]

    else:
        vs0 = self.solver[0,0].getveff()

    #Preallocation

    B = spdiags( 1/n0 , 0, Nelem, Nelem )

    for i in range(self.solver.shape[0]):
        for j in range(self.solver.shape[1]):
            self.solver[i,j].hamiltonian()
            self.solver[i,j].eo = -20

    self.vs = vs0

    print(f"                                         First-order\n")
    print(f"Iteration        f(x)                     optimality\n")

    optimality = 1
    tol        = self.optInv.tol_invert
    maxiter    = self.optInv.max_iter_invert
    iter       = 1

    while optimality > tol and iter < maxiter:
        #Transfer new potentials to solver objects and calculate new densities

        self.solver[i,j].setveff(self.vs)
        self.solver[i,j].calc_orbitals()

        n = np.zeros((self.grid.Nelem, pol))  
        for j in range(self.solver.shape[1]):
            for i in range(self.solver.shape[0]):     
                #Calculate Orbital's densities
                self.solver[i,j].calc_density()
                #Add up orbital's densities
                n[:,j] += np.squeeze(self.solver[i,j].n)

        n = n[:,0]

        grad = B * (n-n0)

        self.vs += 0.1 * grad
        optimality = np.max(np.abs(grad))
        fx         = np.sum(np.abs(grad))
        print(f" {iter}          {fx}           {optimality} ")
        iter += 1 

    max =  0.0 
    for i in range(self.solver.shape[0]):
        for j in range(self.solver.shape[1]):
            self.us = np.max(self.solver[i,j].eig) if np.max(self.solver[i,j].eig) > max else max

    Eks = 0.0

    if tilde is not None:
        self.solver[i,j].calc_energy()

        n = np.zeros((self.grid.Nelem, pol))  
        Eks = 0.0
        for j in range(self.solver.shape[1]):
            for i in range(self.solver.shape[0]):     
                #Calculate Orbital's densities
                self.solver[i,j].calc_density()
                #Add up orbital's densities
                n[:,j] += np.squeeze(self.solver[i,j].n)

        # print("shpae of n", n.shape)
        # print("shape of vs", self.vs[:,None].shape)

        # print("multiply them", n * self.vs[:,None])

        Eks += np.sum(self.solver[i,j].eks)

        self.Ts = Eks - self.grid.integrate(np.sum( n * self.vs[:,None], axis=1 ))



        return True, True
