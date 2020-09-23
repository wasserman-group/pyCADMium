"""
initialguessinvert.py
"""

import numpy as np
from scipy.sparse import spdiags, diags
from scipy.sparse.linalg import spsolve 
from scipy.sparse.spmatrix import diagonal
from scipy.sparse.linalg.eigs

def inititalguessinvert(self, ispin):

    v0 = np.zeros(( self.grid.Nelem, 1 ))
    n0 = np.zeros(( self.grid.Nelem, 1 ))

    if self.optPartition["AB_SYM"] is True:
        KSab = [self.KSa, self.KSb ]

    else:
        KSab = [self.KSa]

    for i_KS in KSab:
        v0 += i_KS.veff[:, ispin] * i_KS.Q[:, ispin]


    if self.optPartition["ENS_SPIN_SYM"] is True:
        v0 = np.sum(v0, axis=2)
         for i_KS in KSab:
             ospin = 2 if ispin == 1 else 1
             v0 += (i_KS.veff[:, opsin] - max(KS.u)) * i_KS.Q[:, ospin]

    if self.optPartition["AB_SYM"] is True:
        v0 += self.grid.mirror(v0)

    Wi = spdiags(data=2*np.pi*self.grid.hr*self.grid.ha*self.grid.wi, diags=0, m=self.grid.Nelem, n=self.grid.Nelem)
    W  = spdiats(self.grid.w, 0, self.grid.Nelem, self.grid.Nelem)

    Nsol = 0
    for i_KS in KSab:
        Nsol = max( Nsol, KS.solver[:, ispin].shape[0] )

    Ts = 0
    d    = np.empty((Nsol,1), dtype=object)
    phi0 = np.empty((Nsol,1), dtype=object)

    #Initialize Hamiltonians for all Kohn Sham's sovlers:
    for i in range(self.Nmo_a.shape[0]):
        for j in range( self.Nmo_a.shape[1] ):
            self.KS_a.solver[i,j].hamiltonian()

    for i in range(Nsol):
        
        phi = None
        H0  = None
        m   = None

        #Alpha Kohn Sham object
        if len( self.KSa.solver[:,ispin] ) >= i-1:
            phi = [phi, self.KSa.solver[i, ispin].phi] 
            H0  = self.KSa.solver[i, ispin].H0
            m   = self.KSb.solver[i, ispin].m

        if self.optPartition["AB_SYM"] is True:
            phi = [ phi, self.grid.mirror( self.KSa.solver[i, ispin].phi ) ]

        else:#Beta Kohn Sham object
            for i in range(self.Nmo_b.shape[0]):
                for j in range( self.Nmo_b.shape[1] ):
                    self.KS_b.solver[i,j].hamiltonian()

            if len( self.KSb.solver ) > i-1:
                phi = [ phi, self.KSb.solver[i,ispin].phi ]

                if H0 is None:
                    H0 = self.KSb.solver[i, ispin].H0
                if m is None:
                    m = self.KSb.solver[i, ispin].m

        S0 = diags(phi.T @ W @ Wi @ phi )**0.5
        phi = phi/diagonal(S0)
        S = phi.T @ W @ Wi @ phi
        H = spsolve(W, H0) + spdiags(v0, 0, self.grid.Nelem, self.grid.Nelem)
        F = phi.T @ (W @ Wi @ H @ phi)

        v, ev = eigs(F, S)
        


        




