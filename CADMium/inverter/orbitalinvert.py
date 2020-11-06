"""
orbitalinvert.py
"""

import numpy as np
import numpy.matlib as matlab

from scipy.sparse import csc_matrix
from scipy.sparse import block_diag as blkdiag
from scipy.sparse import spdiags
from scipy.sparse import eye

import sys
import warnings


def ab_symmetrize(x, Nelem, Nmo, North, self):
    """
    """

    phi   = np.reshape( x[0:Nelem * np.sum(Nmo)], (Nelem, np.sum(Nmo)), order="F")
    evals = np.vstack( (x[ np.array(range(np.sum(Nmo) - 1)) + np.sum(Nmo) * Nelem ], np.array([[0]])) )
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


    for i in range(0, int(np.sum(Nmo, axis=1))):
        if phi[:, i].T @ self.grid.mirror( phi[:, i] ) > 0:
            phi[:, i] += self.grid.mirror( phi[:, i] )
            phi[:, i] /=  np.linalg.norm(  phi[:, i] )

        else:
            phi[:, i] -= self.grid.mirror( phi[:, i] )
            phi[:, i] /= np.linalg.norm( phi[:, i] )

    # x((1:Nelem)+sum(Nmo)*Nelem+sum(Nmo)-1+North) = v;

    x[np.array(range(Nelem * np.sum(Nmo)))] = phi.flatten("F")[:, None]
    v = x[ np.array(range(Nelem)) + np.sum(Nmo) * Nelem + np.sum(Nmo) -1 + North ]
    v = 0.5 * ( v + self.grid.mirror(v))
    x[ np.array(range(Nelem)) + np.sum(Nmo) * Nelem + np.sum(Nmo) -1 + North ] = v
    x = x.real
    
    return x

def normalize(x, Nmo, Nelem, WWi, occ):
    """
    """

    phi = np.reshape( x[ np.array(range(Nelem*np.sum(Nmo))) ], (Nelem, np.sum(Nmo)) , order="F")
    S = np.sum( WWi @ np.abs( phi )**2 , axis=0)
    phi = phi @ np.diag( (occ/S)**0.5 )
    x[0:Nelem * np.sum(Nmo)] = phi.flatten("F")[:, None]


    return x

def kmatrix(x, i_orth, Nelem, Nmo, WWi, H, n0, occ, B2i, B2j, B3i, B3j):

    North = len(i_orth)

    # print("x from x matrix\n", x)
    phi      = np.reshape( x[ np.array(range( Nelem * np.sum(Nmo)))], 
                            (Nelem, np.sum(Nmo)), order="F")
    v        = x[ np.array(range(Nelem)) + np.sum(Nmo) * Nelem + np.sum(Nmo) - 1 + North ]
    evals    = np.vstack(( x[ np.array(range(np.sum(Nmo) - 1)) + np.sum(Nmo) * Nelem ], 0 ))
    if North == 0:
        orthvals = np.zeros((1,0))
    else:
        warnings.warn("Please make sure orthvals are correct")
        orthvals = x[ np.array(range(North)) + np.sum(Nmo) * Nelem + np.sum(Nmo) - 1 ]

    # #BSXFUN
    bsxfun = v - np.ones((v.shape[0], 1)) * evals.T
    vse = WWi @ bsxfun

    Hjac = H + spdiags( vse.flatten('F'), 0 ,np.sum(Nmo) * Nelem, np.sum(Nmo) * Nelem  )
    Ocon = np.zeros((North, 1))

    if North == 0:
        B4 = np.zeros((np.sum(Nmo) * Nelem, 1))
    else:
        B4 = np.zeros(( np.sum(Nmo) * Nelem, North ))
    
    for i in range(North):
        print("Warning North iteration may be *very* wrong")
        Ocon[i] = np.sum( WWi @ phi[:, iorht[i,0]] * phi[:, iorth[i, 1]] )
        ind = np.ravel_multi_index( ( range(0, Nelem) + (i_orth[i, 0]-1)*Nelem ,    
                                      range(0, Nelem) + (i_orth[i, 1]-1)*Nelem)
                                    [np.sum(Nmo) * Nelem, np.sum(Nmo) * Nelem], 
                                    order="F")
                    
        Hjac[ind] = spdiags(WWi) @ orthvals[i]
        ind = np.ravel_multi_index( ( range(0, Nelem) + (i_orth[i, 1]-1)*Nelem ,    
                                      range(0, Nelem) + (i_orth[i, 0]-1)*Nelem)
                                    [np.sum(Nmo) * Nelem, np.sum(Nmo) * Nelem], 
                                    order="F")
        Hjac[ind] = spdiags(WWi) @ orthvals[i]


    KSeq = Hjac @ x[ np.array(range(np.sum(Nmo) * Nelem)) ]

    n = np.sum( np.abs(phi)**2, axis=1)
    S = np.sum( WWi @ np.abs(phi)**2, axis=0)

    ncon = WWi @ ((n-n0) / 2)[:, None]
    Ncon = (- (S - occ) / 2)[:, None]

    # print("Shape of KS", KSeq.shape)
    # print("Ncond", Ncon[:-1].shape)
    # print("ocon", Ocon.shape)
    # print("ncon shape", ncon.shape)

    eqn = np.vstack(( KSeq, Ncon[:-1], Ocon, ncon  ))

    B2v = np.reshape(WWi @ phi, (np.sum(Nmo) * Nelem, 1), order="F")
    B3v = np.reshape(-WWi @ phi[:,:np.sum(Nmo)-1], ((np.sum(Nmo)-1)*Nelem, 1), order="F")

    B2 = csc_matrix( (B2v[:, 0], (B2i[:,0], B2j[:,0])), shape=(np.sum(Nmo) * Nelem, Nelem        ))
    B3 = csc_matrix( (B3v[:, 0], (B3i[:,0], B3j[:,0])), shape=(np.sum(Nmo) * Nelem, np.sum(Nmo)-1))

    B2   = B2.toarray()
    B3   = B3.toarray()
    Hjac = Hjac.toarray()

    if North == 0:
        jac_a = np.concatenate(( Hjac, B3, B2), axis=1)
        jac_b = np.concatenate(( B3.T, B2.T ), axis=0)
    else:
        jac_a = np.concatenate(( Hjac, B3, B4, B2), axis=1)
        jac_b = np.concatenate(( B3.T, B4.T, B2.T ), axis=0)

    rest_matrix = np.zeros((Nelem + np.sum(Nmo) - 1 + North, Nelem + np.sum(Nmo) - 1 + North))
    jac_b = np.concatenate( (jac_b, rest_matrix) , axis=1)
    jac = np.vstack((jac_a, jac_b))

    return jac, eqn


class inversion_info:
    pass

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
    output = {"iterations" : None, "firstorderopt" : None }

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

    if self.vs is None:
        self.vs = np.zeros((self.grid.Nelem, self.solver.shape[1]))

    if self.us is None:
        self.us = np.zeros((self.solver.shape[1]))

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
    for it in range(0, int(np.sum(Nmo))-1 ):
        B3j[ np.array(range(0,Nelem)) + (it) * Nelem ] = it


    #Settings for forcing orthogonality between orbitals
    i_orth = np.zeros( (0, 1) )
    North = len( i_orth )

    #Initial Guess
    if isolver[0].phi is None:
        X = np.vstack((  phi0[:np.sum(Nmo) * Nelem], 
                         (e0[:-2] - e0[np.sum(Nmo)-1])[:, None], 
                         np.zeros((North, 1)),
                         vs0 - e0[np.sum(Nmo)-1]
                        ))


    else:
        
        phi   = (isolver[0].phi)
        evals = (isolver[0].eig)

        for i in range(1,Nsol):
            warnings.warn("WARNING, Double Check concatenation")
            phi   = np.hstack( (phi, (isolver[i].phi) ))
            evals = np.hstack( (evals, isolver[i].eig) )

        # print("phi", phi)
        # print("eig", evals)

        # print("phi", phi.flatten("F").shape)
        # print("evals", (evals[0:-2] - evals[-1]).shape)
        # print("Zeros", np.zeros((North, 1)).shape)
        # print("Pst", (isolver[0].veff - evals[-1]).shape)

        # print("evals", evals)

        X = np.vstack(( phi.flatten("F")[:, None], 
                        (evals[0:-1] - evals[-1]),
                        np.zeros((North, 1)),
                        (isolver[0].veff - evals[-1])[:,None]
                      ))

        # print("x shape", X.shape)


    if self.optInversion["AB_SYM"] is True:
        X = ab_symmetrize(X, Nelem, Nmo, North, self)

    # print("First layer X\n", X)

    X = normalize(X, Nmo, Nelem, WWi, occ)

    # print("Second Layer X\n", X)

    if self.optInversion["DISP"] is True:
        print('iter      res_ks        res_ncon         res_Ncon    res_linsolve  iter_linsolve\n');
        print('--------------------------------------------------------------------------------\n');

    iter        = 0
    res         = 1
    finished    = False
    ForceUpdate = False

    if self.optInversion["AVOIDLOOP"] is True:
        old_res_ks     = 0
        old_res_ncon   = 0
        older_res_ks   = 0
        older_res_ncon = 0

    while finished is not True:
        iter += 1

        ################
        jac, eqn = kmatrix(X, i_orth, Nelem, Nmo, WWi, H, n0, occ, B2i, B2j, B3i, B3j)

        if self.optInversion["AVOIDLOOP"] is True and iter > 1:
            older_res_ks   = old_res_ks
            older_res_ncon = old_res_ncon
            old_res_ks     = res_ks
            old_res_ncon   = res_ncon

        # print("Ecuation\n", eqn.T)



        res_ks   = np.max(np.abs( eqn[ np.array(range(Nelem * np.sum(Nmo))) ]))
        res_ncon = np.max(np.abs(
                    eqn[ np.array(range( Nelem )) + np.sum(Nmo) * Nelem + np.sum(Nmo) - 1 + North  ]
                    ))
        res_Ncon = np.max(np.abs( eqn[ np.array(range(np.sum(Nmo)-1)) +Nelem * np.sum(Nmo) ] ))
        res = max(res_ks, res_ncon)


        if self.optInversion["AVOIDLOOP"] is True and iter > 2:
            if res      < 1e-6 and \
               res_ks   > old_res_ks     * 0.1 and \
               res_ks   < older_res_ks   * 10. and \
               res_ncon > older_res_ncon * 0.1 and \
               res_ncon < older_res_ncon * 10.:

                ForceUpdate = True
            else:
                ForceUpdate = False

        if self.optInversion["DISP"]:
             print(f"{iter} {res_ks} {res_ncon} {res_Ncon}")

        if res >1e2:
            print("ERROR too large")
            finished = True

        # print(f"{res_ncon} >? {res_ks * self.optInversion['ResFactor']}")
        # print(f"{res_ks}   <? {self.optInversion['TolInvert'] * 1e2}")
        # # print(f"{ForceUpdate}")

        # print(f"{res_ncon > (res_ks * self.optInversion['ResFactor'])}, \
        #         {res_ks   < (self.optInversion['TolInvert'] * 1e2)},    \
        #         {ForceUpdate}")


        if res < self.optInversion["TolInvert"]:
            finished = True

        elif iter >= self.optInversion["MaxIterInvert"]:
            finished = True
            warnings.warn('\n Convergence not reached at maximum iteration')

        else:
            #Convergence restrictions have been met
            if (res_ncon > (res_ks * self.optInversion["ResFactor"]) or \
                res_ks   < (self.optInversion["TolInvert"] * 1e2) or    \
               ForceUpdate) == True:

                #Use Iterative?
                if self.optInversion["USE_ITERATIVE"] is True:
                    raise ValueError("USE ITERATIVE Flag Not Avaliable")
                    # tol = self.optInversion["Tolinsolve"]
                    # warnings.warn("Nearly Singular Matrix")

                    # A = jac[:np.sum(Nmo) * Nelem, :np.sum(Nmo) * Nelem]
                    # B = jac[:np.sum(Nmo) * Nelem, np.sum(Nmo) * Nelem + 1:]

                    # M1 = np.concatenate((A, B), axis=1)
                    # M2 = np.concatenate((B.T, eye(B.shape[1])), axis=1)
                    # M = np.concatenate( (M1, M2), axis=0 )
                #Or not use iterative
                else:
                    #If use iterative is False
                    dX = - np.linalg.solve(jac, eqn)

                #Add dX to X
                dv = dX[ np.array(range( Nelem )) + np.sum(Nmo) * Nelem + np.sum(Nmo) -1 + North ]
                X += dX

            #If conditions have not been met, recalculate X
            else:
                self.vs[:, ispin] = X[np.array(range(Nelem)) + np.sum(Nmo) * Nelem + np.sum(Nmo) - 1 + North, 0]
                evals = np.concatenate( ( X[ np.array(range(np.sum(Nmo) - 1)) + np.sum(Nmo) * Nelem ], np.array(([[0]])) ), axis=0 )

                for it in range(Nsol):
                    isolver[it].phi = np.reshape( X[ np.array(range(Nmo[it] * Nelem)) + np.sum(Nmo[0:it-1]) * Nelem] , (Nelem, Nmo[it]) , order="F")
                    isolver[it].eig = evals[:Nmo[it] + np.sum(Nmo[:it-1])]

                for i in isolver:
                    i.setveff(self.vs[:, ispin])
                    i.calc_orbitals()

                phi = None
                evals = None

                for i in isolver:
                    if phi is None:
                        phi = i.phi
                    else:
                        phi   = np.concatenate( (phi, i.phi),   axis=1)

                    if evals is None:
                        evals = i.eig
                    else:
                        evals = np.concatenate( (evals, i.eig), axis=1)


                X = np.concatenate((phi.flatten("F")[:,None], 
                                    (evals[:-1] - evals[-1])[:, None], 
                                    np.zeros((North, 1)), 
                                    (self.vs[:, ispin]-evals[-1])[:, None]),
                                    axis=0 )
                        
            if self.optInversion["AB_SYM"] is True:
                X = ab_symmetrize(X, Nelem, Nmo, North, self)

            X = normalize(X, Nmo, Nelem, WWi, occ)

        if self.optInversion["DISP"] is True:
            pass


    self.vs[:, ispin] = X[ np.array(range(Nelem)) + np.sum(Nmo) * Nelem + np.sum(Nmo) - 1 + North, 0]
    evals = np.concatenate( ( X[ np.array(range(np.sum(Nmo) - 1)) + np.sum(Nmo) * Nelem ], np.array(([[0]])) ), axis=0 )

    for it in range(Nsol):
        isolver[it].phi = np.reshape( X[ np.array(range(Nmo[it] * Nelem)) + np.sum(Nmo[0:it-1]) * Nelem] , (Nelem, Nmo[it]) , order="F")
        isolver[it].eig = evals[:Nmo[it] + np.sum(Nmo[:it-1])]

    for it in range(Nsol):
        self.us[ispin] = isolver[it].get_homo()

        isolver[it].calc_density()
        isolver[it].setveff(self.vs[:, ispin])

    if res > self.optInversion["TolInvert"]:
        flag = False
    else:
        flag = True


    inv_info = inversion_info
    inv_info.nfev       = iter 
    inv_info.optimality = res
    output = np.empty((1,1), dtype=object)
    output[0,0] = inv_info

    return flag, output







