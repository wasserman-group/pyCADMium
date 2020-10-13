"""
calc_response.py
"""

import numpy as np
from scipy.sparse import spdiags
from scipy.sparse import vstack as csc_vstack
from scipy.sparse import hstack as csc_hstack
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import splu

def calc_response(self):
    """
    Calculate the response function
    """

    Nelem = self.grid.Nelem

    #Inverse volume element
    W = spdiags(self.grid.w, 0, Nelem, Nelem)

    #Build effective potential operator
    Veff = spdiags(W @ self.veff, 0, Nelem, Nelem)

    self.chi = np.zeros((Nelem))    
    H = self.H0 + Veff

    #Figure out the number of occupied norbitals: Nocc
    if self.m == 0:
        if self.pol == 1:
            Nocc = np.floor(self.N/2)
            nu = self.N / 2 - Nocc
        else:
            Nocc = np.floor(self.N)
            nu = self.N - Nocc

    else:
        #m>0 orbitals hold twice as many electrons due to +-m symmetry
        if self.pol == 1:
            Nocc = np.floor(self.N / 4)
            nu = self.N / 4 - Nocc
        else:
            Nocc = np.floor(self.N/2)
            nu = self.N / 2 - Nocc 

    #Calculate response from fully occupied orbitals
    for j in range(int(Nocc)):

        #Select an orbital and normalize
        phi_j = self.phi[:, j] / self.grid.integrate( self.phi[:, j]**2)**0.5

        #Construct right hand side
        PHI = spdiags(phi_j, 0, Nelem, Nelem)
        rhs = csc_vstack( ( W @ PHI, np.zeros((1, Nelem))))
        rhs = csc_matrix(rhs)
        G = H - W * self.eig[j]
    
        A_inside_1 = csc_hstack(( -G, (W @ phi_j)[:,None] ))
        A_inside_2 = csc_hstack(( (W * phi_j).T, np.array([[0]]) ))
        A = csc_vstack(( A_inside_1, A_inside_2 ))
        A = csc_matrix(A)

        x = spsolve(A, rhs)

        #Add orbital density response into total resposne
        self.chi += 2 * PHI @ x[0:Nelem, 0:Nelem]

    #If we are doing a fractional orbitals and are non-integer
    if self.FRACTIONAL is True and nu != 0:
        print("Warning, you may need to confirm calc_response in CADMium")
        j = Nocc + 1

        #Select the orbital and normalize
        phi_j = self.phi[:, j] / self.grid.integrate( self.phi[:, j]**2 )**0.5

        #Construct the right hand side
        PHI = spdiags(phi_j, 0, Nelem, Nelem)
        rhs = csc_vstack( ( W @ PHI, np.zeros((1, Nelem))))
        rhs = csc_matrix(rhs)

        #Left hand side
        G = H - W * self.eig[j]
        A_inside_1 = csc_hstack(( -G, (W @ phi_j)[:,None] ))
        A_inside_2 = csc_hstack(( (W * phi_j).T, np.array([[0]]) ))
        A = csc_vstack(( A_inside_1, A_inside_2 ))
        A = csc_matrix(A)

        lu = splu(A)
        #Take matrox L and U from decomposition
        L = lu.L.A
        U = lu.U.A

        #Solve linear equation
        x = spsolve( U , spsolve(L, rhs))

        #Add orbital density response into total response
        self.chi += self.chi + nu * 2 * PHI * x[0:Nelem, 0:Nelem]

    #Scale responses accordingly
    if self.pol == 1:
        self.chi = 2 * self.chi


