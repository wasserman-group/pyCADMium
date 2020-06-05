"""

Operators.py
Construct operators

"""

import numpy as np
from scipy.sparse import eye, csc_matrix, kron, spdiags
#from .mirror import mirror

def operators(self):

    #Scale factors 
    self.h1 = self.a * (np.sinh(self.Xr)**2 + np.sin(self.Xa)**2)**0.5
    self.h2 = self.a * (np.sinh(self.Xr)**2 + np.sin(self.Xa)**2)**0.5
    self.h3 = self.a * np.sinh(self.Xr) * np.sin(self.Xa)
    
    self.h1 = 0.5 * (self.h1 + self.mirror(self.h1))
    self.h2 = 0.5 * (self.h2 + self.mirror(self.h2))
    self.h3 = 0.5 * (self.h3 + self.mirror(self.h3))

    #Finite difference combined operators
    Da1 = kron(eye(self.Nr), csc_matrix(self.eDa1), format="csc")
    Dr1 = kron(csc_matrix(self.eDr1), eye(self.Na), format="csc")

    self.gradr = spdiags(data = 1.0 / self.h2, diags = 0, m = self.Nelem, n = self.Nelem) * Dr1
    self.grada = spdiags(data = 1.0 / self.h1, diags = 0, m = self.Nelem, n = self.Nelem) * Da1

    #Finite difference combined operators
    Da1 = kron(eye(self.Nr), csc_matrix(self.oDa1), format="csc")
    Dr1 = kron(csc_matrix(self.oDr1), eye(self.Na), format="csc")

    C123 = spdiags(data = 1.0 / (self.h1 * self.h2 * self.h3), diags = 0, m = self.Nelem, n = self.Nelem)
    self.diva = C123 * Da1 * spdiags(data = self.h2 * self.h3, diags = 0, m = self.Nelem, n = self.Nelem)
    self.divr = C123 * Dr1 * spdiags(data = self.h1 * self.h3, diags = 0, m = self.Nelem, n = self.Nelem)

    #Construct 1D grids
    Xa1d = np.reshape(self.Xa, (self.Na, self.Nr), order="F")
    Xr1d = np.reshape(self.Xr, (self.Na, self.Nr), order="F")
    Xa1d = Xa1d[:, 0][None].T
    Xr1d = Xr1d[0, :][None]

    
    #Even Symmetry
    #Angular and radial portions of laplacian
    eLa = self.a * (np.diagflat(np.sin(Xa1d))@self.eDa2 + np.diagflat(np.cos(Xa1d))@self.eDa1)
    eLr = self.a * (np.diagflat(np.sinh(Xr1d)) @ self.eDr2 + np.diagflat(np.cosh(Xr1d)) @ self.eDr1)
    #Odd Symmetry
    #Angular and radial portions of laplacian
    oLa = self.a * (np.diagflat(np.sin(Xa1d))@self.oDa2 + np.diagflat(np.cos(Xa1d))@self.oDa1)
    oLr = self.a * (np.diagflat(np.sinh(Xr1d))@self.oDr2 + np.diagflat(np.cosh(Xr1d))@self.oDr1)


    #Even Symmetry
    #Laplacian
    self.elap = (kron(csc_matrix(np.diagflat(np.sinh(Xr1d))), csc_matrix(eLa), format="csc") 
                + kron(csc_matrix(eLr), csc_matrix(np.diagflat(np.sin(Xa1d))), format="csc"))

    #Odd Symmetry
    #Laplacian
    self.olap = (kron(csc_matrix(np.diagflat(np.sinh(Xr1d))), csc_matrix(oLa), format="csc") 
               + kron(csc_matrix(oLr), csc_matrix(np.diagflat(np.sin(Xa1d))), format="csc"))



    #Construct Laplacian operator for boundary values
    bLr = self.a * (  np.diagflat(np.sinh(Xr1d[0][-self.bcN:])) @ self.bc2  
                    + np.diagflat(np.cosh(Xr1d[0][-self.bcN:])) @ self.bc1)

    self.blap = kron(csc_matrix(bLr), csc_matrix(np.diagflat(np.sin(Xa1d))), format="csc")

