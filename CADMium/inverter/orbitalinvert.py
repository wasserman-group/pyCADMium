"""
orbitalinvert.py
"""

import numpy as np
import numpy.matlib as matlab

from scipy.sparse import block_diag as blkdiag
from scipy.sparse import spdiags

import sys

def ab_symmetrize(x, Nelem, Nmo, North, self):
    """
    """

    phi   = np.reshape( x[0:Nelem * np.sum(Nmo)], (Nelem, np.sum(Nmo)), order="F")
    evals = np.vstack( (x[0 : (np.sum(Nmo) - 2) + np.sum(Nmo) * Nelem], np.array([[0]])) )

    Nmo = Nmo[:, None]

    #Check for degenerate and nearly degenerate orbitals
    for i in range(0, int(np.sum(Nmo, axis=1)-2)):
        if np.abs( evals[i] - evals[i+1] ) < 1e-6:

            even =  phi[:, i]   + self.grid.mirror( phi[:,i]   ) \
                  + phi[:, i+1] + self.grid.mirror( phi[:,i+1] )
            odd  =  phi[:, i]   - self.grid.mirror( phi[:,i]   ) \
                  + phi[:, i+1] - self.grid.mirror( phi[:,i+1] )

            phi[:, i]   = even / np.linalg.norm(even)
            phi[:, i+1] = odd  / np.linalg.norm(odd)

    for i in range(0, int(np.sum(Nmo, axis=1) - 1)):
        if phi[:, i].T @ self.grid.mirror( phi[:, i] ) > 0:
            phi[:, i] = phi[:, i] + self.grid.mirror( phi[:,i] )
            phi[:, i] = phi[:, i] / np.linalg.norm( phi[:, i] )

        else:
            phi[:, i] = phi[:, i] - self.grid.mirror( phi[:, i] )
            phi[:, i] = phi[:, i] / np.linalg.norm( phi[:, i] )

    x[0:Nelem * np.sum(Nmo)] = phi.flatten("F")[:, None]
    v = x[ np.array(range(0, Nelem)) + np.sum(Nmo) * Nelem + np.sum(Nmo) - 1 + North ]
    v = 0.5 * ( v + self.grid.mirror(v))
    x[ np.array(range(0, Nelem)) + np.sum(Nmo) * np.sum(Nmo) -1 + North] = v
    x = x.real

    return x

def normalize(x, Nmo, Nelem, WWi, occ):
    """
    """

    phi = np.reshape( x[0:Nelem*np.sum(Nmo)], (Nelem, np.sum(Nmo)) )
    S = np.sum( WWi @ np.abs( phi )**2 )
    phi = phi @ np.diag( (occ/S)**0.5 )
    x[0:Nelem * np.sum(Nmo)] = phi.flatten("F")[:, None]

    return x

def orbitalinvert(self, n0, vs0, phi0, e0, ispin):
    """
    Orbital Invert

    Jacobian:
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %          %      %    %
    %          %      %    %
    %   Hjac   %  B3  % B2 %
    %          %      %    %
    %          %      %    %
    %%%%%%%%%%%%%%%%%%%%%%%%
    %          %      %    %
    %    B3'   %  0   % 0  %
    %          %      %    %
    %%%%%%%%%%%%%%%%%%%%%%%%
    %    B2'   %  0   % 0  %
    %%%%%%%%%%%%%%%%%%%%%%%%

    """

    Nelem = self.grid.Nelem
    Wi = spdiags( 2*np.pi * self.grid.hr * self.grid.ha * self.grid.wi, 0, Nelem, Nelem)
    W  = spdiags( self.grid.w, 0, Nelem, Nelem )
    WWi = W @ Wi

    isolver = self.solver[:, ispin]
    Nsol = isolver.shape[0]
    #Nmo, N, m, = [], [], []
    Nmo, N, m = np.empty((Nsol), dtype=int), np.empty((Nsol), dtype=int), np.empty((Nsol), dtype=int)
    for i, j in enumerate(isolver):
        j.hamiltonian()
        Nmo[i] = j.Nmo
        N[i] = j.N
        m[i] = j.m
    Nsol = isolver.shape[0]

    #Transforming cell array from matlab to numpy array
    C = np.empty((Nsol), dtype=object)
    D = np.empty((Nsol), dtype=object)

    for it in range(Nsol):
        C[it] = matlab.repmat( Wi @ isolver[it].H0, Nmo[it], 1 )

        if m[it] == 0:            
            D[it] = 2 * np.ones( (Nmo[it], 1) )

            if isolver[it].pol == 1:
                nu = N[it] / 2 - np.floor(N[it] / 2)
            else:
                nu = N[it] - np.floor(N[it])
            
            if nu != 0:
                D[it][Nmo[it]] = D[it][Nmo[it]] * nu
        
        else:
            D[it] = 4 *  np.ones((Nmo[it,0], 1))
            if isolver[it].pol == 1:
                nu = N[it] / 4 - np.floor(N[it] / 4)
            else:
                nu = N[it] / 2 - np.floor(N[it] / 2)

            if nu != 0:
                D[it][Nmo[it]] = D[it][Nmo[it]] * nu

    C = np.squeeze(np.vstack((C[:])))
    H = blkdiag(  C[:] )
    occ = np.squeeze(np.vstack(D[:]))

    B2i = np.array( range( 0, np.sum(Nmo) * Nelem ) )[:, None]
    B2j = matlab.repmat( np.array(range(0,Nelem))[:,None], np.sum(Nmo), 1 )
    B3i = np.array( range( 0,(np.sum(Nmo) - 1) * Nelem ) )[:, None]
    B3j = np.zeros( ((np.sum(Nmo) - 1) * Nelem, 1) )

    for it in range( 1, int(np.sum(Nmo) + 0) ):
        B3j[ np.linspace(0,Nelem-1, Nelem, dtype=int) + (it - 1) * Nelem ] = it


    #Settings for forcing orthogonality between orbitals
    i_orth = []
    North = len( i_orth )

    #Initial Guess
    if isolver[0].phi is None:
        X = np.vstack((  phi0[np.array(range(0,np.sum(Nmo) * Nelem) ,dtype=int)][:, None], 
                         (e0[0:np.sum(Nmo)-1] - e0[np.sum(Nmo)])[:, None], 
                         np.zeros((North, 1)),
                         vs0 - e0[np.sum(Nmo)]
                        ))

    else:
        print("WARNING, Concatenation will probably be wrong")
        phi   = isolver[0].phi[:,0]
        evals = isolver[0].eig[0]
        for i in range(1,Nsol):
            phi   = np.hstack( (phi, isolver[i].phi[:,i]) )
            evals = np.hstack( (evals, isolver[i].eig[i]) )
        X = np.vstack(( phi, 
                        evals[0:-2] - evals[-1],
                        np.zeros((North, 1)),
                        isolve[0].veff - evals[-1]
                      ))

    if self.optInversion["AB_SYM"] is True:
        X = ab_symmetrize(X, Nelem, Nmo, North, self)

    X = normalize(X, Nmo, Nelem, WWi, occ)

    print(X)
        

    sys.exit()







