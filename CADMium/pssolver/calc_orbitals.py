"""
calc_orbitals.py
"""

from scipy.sparse import spdiags
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import eigs
from scipy.sparse.linalg import spsolve
from numpy.linalg import norm
from numpy import spacing
import numpy as np

import sys

def calc_orbitals(self, solver_id, return_dict):
    """
    calculate molecular orbitals and eigenvalues
    """

    assert len(self.veff) != 0, "Veff is not set"

    Nelem = self.grid.Nelem

    if self.Nmo != 0:

        #Inverse volume element
        W = spdiags(data=self.grid.w, diags=0, m=Nelem, n=Nelem)
        W = csc_matrix(W)

        #Build effective potential operator
        Veff = spdiags(data= W @ self.veff, diags=0, m=Nelem, n=Nelem)

        if self.H0 is None:
            self.hamiltionian()

        #Construct Hamiltonian
        H = self.H0 + Veff

        #Solve eigenvalue problem
        eig, phi = eigs(spsolve(W, H), k=self.Nmo, sigma=self.e0, v0=self.opt["v0"])
        eig = eig.real
        phi = phi.real
        e0 = self.e0

        #Ordering Orbitals. 
        indx = eig.argsort()
        self.eig = eig[indx]
        self.phi = phi[:,indx] 

        while np.isnan(self.phi).all() != np.zeros_like(self.phi).all():
            e0 = e0 - 0.1
            self.eig, self.phi = eigs(spsolve(W, H), k=self.Nmo, sigma=e0, v0=self.opt["v0"])
            eig = self.eig.real
            phi = self.phi.real

        #Check for degenerate and nearly degenerate orbitals
        for i in range(self.Nmo-1):
            for j in range(i+1, self.Nmo):
                if np.abs(self.eig[i]-self.eig[j]) < 1e-9:
                    even = self.phi[:, i] + self.grid.mirror(self.phi[:,i]) + self.phi[:,j] + self.grid.mirror(self.phi[:,j])
                    odd  = self.phi[:, i] - self.grid.mirror(self.phi[:,i]) + self.phi[:,j] - self.grid.mirror(self.phi[:,j])
                    self.phi[:, i] = even/norm(even)
                    self.phi[:, j] = odd/norm(odd)

        if self.optSolver.sym is True:
            for i in range(self.Nmo):
                if self.phi[:,i].T @ self.grid.mirror(self.phi[:,i]) > 0:

                    self.phi[:,i] = self.phi[:, i] + self.grid.mirror(self.phi[:,i])
                    self.phi[:,i] = self.phi[:, i] / norm(self.phi[:, i])

                else:

                    self.phi[:, i] = self.phi[:, i] - self.grid.mirror(self.phi[:,i])

    else:
        self.eig = -1.0 / np.finfo(float).eps
    if return_dict is not None:
        return_dict[solver_id] = [self.eig, self.phi] 