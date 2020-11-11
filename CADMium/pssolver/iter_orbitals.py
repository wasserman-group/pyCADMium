"""
iter_orbitals.py
"""

import numpy as np
from scipy.sparse import spdiags
from scipy.sparse import csc_matrix
from scipy.sparse import hstack
from scipy.sparse import vstack
from scipy.sparse.linalg import spilu
from scipy.sparse.linalg import bicgstab
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import LinearOperator


def iter_orbitals(self, solver_id, return_dict):
    """
    Update molecular orbitals and eigenvalues iteratively
    """
    if self.Nmo != 0:

        assert len(self.veff) != 0, "Veff is not set"
        if self.phi is None or self.eig is None:
            self.calc_orbitals()

        Nelem = self.grid.Nelem
        W = spdiags(self.grid.w, 0, Nelem, Nelem)

        #Construct Hamiltonian
        Hks = self.H0 + spdiags(self.grid.w * self.veff, 0, Nelem, Nelem)
        resvec = Hks @ self.phi - W @ (self.phi @ np.diag(self.eig))
        res = np.amax(np.abs(resvec), axis=0)

        for j in range(self.Nmo):
            if res[j] > self.TOL_ORBITAL:
                uno = Hks - W * self.eig[j]
                dos =  (-W @ self.phi[:,j])[None].T
                tres = -(W @ self.phi[:,j]).T
                cuatro = np.array([-1e-12])

                A = hstack((uno,dos))
                B = hstack((tres, cuatro))
                C = vstack((A,B))

                rhs = np.hstack((resvec[:, j], [0]))
                C = csc_matrix(C)

                if self.ITERLINSOLVE is True:
                    ILU = spilu(C)
                    approx_sol = LinearOperator((C.shape[0], C.shape[1]), ILU.solve)
                    x = bicgstab(C, rhs, tol=1e-15, M=approx_sol)[0]
                
                else:
                    x = spsolve(C, rhs)
                
            self.phi[:, j] -= x[:Nelem]
            self.eig[j]    -= x[-1]
            self.normalize_orbitals()
        
    else:
        self.eig = -1 / np.spacing(1)

    return_dict[solver_id] = [self.eig, self.phi] 
