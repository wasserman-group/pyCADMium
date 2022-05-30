"""
initialguessinvert.py
"""

import numpy as np
from scipy.sparse import spdiags, csc_matrix
from scipy.sparse.linalg import spsolve 
from scipy.linalg import eig, eigh
# from scipy.sparse.linalg import eigs

def initialguessinvert(self, ispin=0):

    v0 = np.zeros(( self.grid.Nelem, 1 ))
    # n0 = np.zeros(( self.grid.Nelem, 1 ))

    if self.optPartition.ab_sym is not True:
        if not self.ens:
            KSab = [self.KSa, self.KSb]
        else:
            KSab = [self.KSa, self.KSA, self.KSb, self.KSB]
    else:
        if not self.ens:
            KSab = [self.KSa]
        else:
            KSab = [self.KSa, self.KSA]

    for i_KS in KSab:
        v0 += i_KS.veff[:, ispin][:, None] * i_KS.Q[:, ispin][:, None]

    if self.optPartition.ens_spin_sym:
        v0 = np.sum(v0, axis=1)
        for i_KS in KSab:
            ospin = 1 if ispin == 0 else 0
            v0 += (i_KS.veff[:, ospin] - max([i_KS.u])) * i_KS.Q[:, ospin]
        v0 = v0[:, None]

    if self.optPartition.ab_sym is True:
        v0 += self.grid.mirror(v0) 


    Wi = spdiags(data=2*np.pi*self.grid.hr*self.grid.ha*self.grid.wi, diags=0, m=self.grid.Nelem, n=self.grid.Nelem)
    W  = spdiags(self.grid.w, 0, self.grid.Nelem, self.grid.Nelem)

    Nsol = 0
    for i_KS in KSab:
        Nsol = max( Nsol, i_KS.solver[:, ispin].shape[0] )

    # Ts = 0
    # Stores the final product
    # d    = []
    # phi0 = []
    phi0 = np.empty( (self.grid.Nelem, 0) )
    d    = np.empty( (              0, 0) )


    for i in range(1, Nsol+1):
        phi = np.empty((self.grid.Nelem,0))
        H0  = None
        m   = None

        #--> Alpha Kohn Sham object
        if len( self.KSa.solver[:,ispin]) >= i-1:
            phi = np.append(phi, self.KSa.solver[i-1, ispin].phi, axis = 1 )
            self.KSa.solver[i-1, ispin].hamiltonian()
            H0  = self.KSa.solver[i-1, ispin].H0
            m   = self.KSb.solver[i-1, ispin].m

        if self.ens: # Ensemble?
            if len( self.KSA.solver[:,ispin]) >= i-1:
                phi = np.append(phi, self.KSA.solver[i-1, ispin].phi, axis = 1)
                self.KSA.solver[i-1, ispin].hamiltonian()
                H0  = self.KSA.solver[i-1, ispin].H0
                m   = self.KSA.solver[i-1, ispin].m

        if self.optPartition.ab_sym is True:
            phi = np.append(phi, self.KSa.solver[i-1, ispin].phi, axis = 1 )
            if self.ens:
                phi = np.append(phi, self.KSA.solver[i-1, ispin].phi, axis = 1 )

        #--> Beta Kohn Sham object
        else:
            if len(self.KSb.solver) >= i-1:
                phi = np.append(phi, self.KSb.solver[i-1, ispin].phi, axis = 1 )
                if H0 is None:
                    H0 = self.KSb.solver[i-1, ispin].H0
                if m is None:
                    m = self.KSb.solver[i-1, ispin].m
            
            if self.ens: # Ensemble?
                if len(self.KSB.solver) >= i-1:
                    phi = np.append(phi, self.KSB.solver[i-1, ispin].phi, axis = 1 )
                    if H0 is None:
                        H0 = self.KSB.solver[i-1, ispin].H0
                    if m is None:
                        m = self.KSB.solver[i-1, ispin].m
        
        S0 = np.diag(phi.T @ W @ Wi @ phi) ** 0.5
        phi = phi / S0
        phi = np.nan_to_num(phi, nan=0.0, posinf=0.0, neginf=0.0)
        S = phi.T @ W @ Wi @ phi
        H = spsolve(csc_matrix(W), H0) + csc_matrix(spdiags(v0[:,0], 0, self.grid.Nelem, self.grid.Nelem))
        F = phi.T @ (W @ Wi @ H @ phi)
        ev, v = eig(a=F, b=S)
        indx = np.argsort(ev)
        d = np.append( d, ev[indx].real )
        v = v[:, indx]
        phi = phi @ v


        if self.optPartition.ens_spin_sym:
            if m == 0:
                phi_to_zero = (phi / np.diag(phi.T @ W @ Wi @ phi)**0.5) * (1**0.5)
                phi0 = np.append( phi0, phi_to_zero , axis = 1)
                # phi0.append( phi_to_zero.real )
                # Eks  = np.sum(d)
            else:
                phi_to_zero = (phi / np.diag(phi.T @ W @ Wi @ phi)**0.5) * (2**0.5)
                phi0 = np.append( phi0, phi_to_zero , axis = 1 )
                # phi0.append( phi_to_zero.real )
                # Eks = 2 * np.sum(d)

        else:
            if m == 0:
                phi_to_zero = (phi / np.diag(phi.T @ W @ Wi @ phi)**0.5) * (2**(0.5))
                phi0 = np.append( phi0, phi_to_zero , axis =  1)
                # phi0.append( phi_to_zero.real )
                # Eks = 2 * np.sum(d)
            else: 
                phi_to_zero = (phi / np.diag(phi.T @ W @ Wi @ phi)**0.5) * (4**(0.5))
                phi0 = np.append( phi0, phi_to_zero , axis = 1)
                # phi0.append( phi_to_zero.real )
                # Eks = 4 * np.sum(d)

        # Ts += Eks - self.grid.integrate( np.sum(phi0**2, axis=1) * np.squeeze(v0) )
        # n0 += np.sum( phi0**2, axis=1 )[:, None]

    e0   = np.nan_to_num(    d, nan=0.0, posinf=0.0, neginf=0.0 )
    phi0 = np.nan_to_num( phi0, nan=0.0, posinf=0.0, neginf=0.0 )

    return phi0, e0, v0


        




