"""
orbitalinvert.py
"""

import numpy as np
import numpy.matlib as matlab

from scipy.sparse import block_diag as blkdiag
from scipy.sparse import spdiags

import sys

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
    Nmo, N, m, = [], [], []
    for i in isolver:
        i.hamiltonian()
        Nmo.append(i.Nmo)
        N.append(i.N)
        m.append(i.m)
    Nsol = isolver.shape[0]

    #Transforming cell array from matlab to numpy array
    C = np.empty((Nsol,1), dtype=object)
    D = np.empty((Nsol,1), dtype=object)

    for it in range(Nsol):
        C[it,0] = matlab.repmat( Wi @ isolver[it].H0, Nmo[it], 1 )

        if m[it] == 0:            
            D[it,0] = 2 * np.ones( (Nmo[it], 1) )

            if isolver[it].pol == 1:
                nu = N[it] / 2 - np.floor(N[it] / 2)
            else:
                nu = N[it] - np.floor(N[it])
            
            if nu != 0:
                D[it,0][Nmo[it]] = D[it,0][Nmo[it]] * nu
        
        else:
            D[it,0] = 4 *  np.ones((Nmo[it,0], 1))
            if isolver[it].pol == 1:
                nu = N[it] / 4 - np.floor(N[it] / 4)
            else:
                nu = N[it] / 2 - np.floor(N[it] / 2)

            if nu != 0:
                D[it,0][Nmo[it]] = D[it,0][Nmo[it]] * nu

    print("This is C", C[0,0])

    print(C[:,0])

    C = np.vstack( (C[:]))
    H = blkdiag(  C[:] )
    occ = np.vstack(D[:] )
            

    sys.exit()



