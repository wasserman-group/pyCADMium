"""
orbitalinvert.py
"""

import numpy as np
import numpy.matlib as matlab

from scipy.sparse import csc_matrix
from scipy.sparse import block_diag
from scipy.sparse import spdiags
from scipy.sparse import eye
from scipy.sparse.linalg import qmr
from scipy.sparse.linalg import LinearOperator
#from scipy.sparse.linalg import spilu
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import inv
#from scipy.sparse.linalg import lgmres

import warnings

class inversion_info:
    pass


def ab_symmetrize(x, Nelem, Nmo, North, self):
    """
    """

    #Inner range of phi:
    if np.sum(Nmo) - 1 == 0:
        evals = np.vstack( (np.array([[0]])) )
    else:
        subset = np.array(range(np.sum(Nmo) - 1)) + np.sum(Nmo) * Nelem
        evals = np.vstack( (x[subset], np.array([[0]])) )

    phi   = np.reshape( x[0:Nelem * np.sum(Nmo)], (Nelem, np.sum(Nmo)), order="F")
    #evals = np.vstack( (x[subset], np.array([[0]])) )
    Nmo = Nmo[:, None]

    #Check for degenerate and nearly degenerate orbitals
    for i in range(0, int(np.sum(Nmo, axis=1)-1)):
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
    S = np.sum( WWi @ np.abs( phi )**2 , axis=0)[:, None]
    phi = phi @ np.diag( ((occ/S)**0.5)[:,0] )
    phi = phi.flatten("F")
    x[0:Nelem * np.sum(Nmo)] = phi[:, None]

    return x

def kmatrix(x, i_orth, Nelem, Nmo, WWi, H, n0, occ, B2i, B2j, B3i, B3j):

    North = len(i_orth)

    #Inner range of phi:
    if np.sum(Nmo) - 1 == 0:
        evals    = np.array([[0]])
    else:
        subset = np.array(range(np.sum(Nmo) - 1)) + np.sum(Nmo) * Nelem
        evals    = np.vstack(( x[ subset ], 0 ))

    # print("x from x matrix\n", x)
    phi      = np.reshape( x[ np.array(range( Nelem * np.sum(Nmo)))], 
                            (Nelem, np.sum(Nmo)), order="F")
    
    v        = x[ np.array(range(Nelem)) + np.sum(Nmo) * Nelem + np.sum(Nmo) - 1 + North ]


    if North == 0:
        orthvals = np.zeros((1,0))
    else:
        warnings.warn("Please make sure orthvals are correct")
        orthvals = x[ np.array(range(North)) + np.sum(Nmo) * Nelem + np.sum(Nmo) - 1 ]


    # #BSXFUN
    bsxfun = v - np.ones((v.shape[0], 1)) * evals.T
    #bsxfun = v - np.dot(v, evals)
    vse = WWi @ bsxfun
    Hjac = H + spdiags( vse.flatten('F'), 0 ,np.sum(Nmo) * Nelem, np.sum(Nmo) * Nelem  )
    Ocon = np.zeros((North, 1))


    if North == 0:
        B4 = np.zeros((np.sum(Nmo) * Nelem, 1))
    else:
        B4 = np.zeros(( np.sum(Nmo) * Nelem, North ))
    
    for i in range(North):
        raise ValueError("pyCADMium Error. Exit")
        # print("Warning North iteration may be *very* wrong")
        # Ocon[i] = np.sum( WWi @ phi[:, iorht[i,0]] * phi[:, iorth[i, 1]] )
        # ind = np.ravel_multi_index( ( range(0, Nelem) + (i_orth[i, 0]-1)*Nelem ,
        #                               range(0, Nelem) + (i_orth[i, 1]-1)*Nelem),
        #                             [np.sum(Nmo) * Nelem, np.sum(Nmo) * Nelem], 
        #                             order="F")
                    
        # Hjac[ind] = spdiags(WWi) @ orthvals[i]
        # ind = np.ravel_multi_index( ( range(0, Nelem) + (i_orth[i, 1]-1)*Nelem , 
        #                               range(0, Nelem) + (i_orth[i, 0]-1)*Nelem),
        #                             [np.sum(Nmo) * Nelem, np.sum(Nmo) * Nelem], 
        #                             order="F")


    KSeq = Hjac @ x[ np.array(range(np.sum(Nmo) * Nelem)) ]

    n = np.sum( np.abs(phi)**2, axis=1)
    S = np.sum( WWi @ np.abs(phi)**2, axis=0)
    ncon = WWi @ ((n-n0) / 2)[:, None]
    Ncon = (- (S - occ) / 2)[:,0][:,None]
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

    if self.vs is None and self.us is None:
        self.vs = np.zeros((self.grid.Nelem, self.solver.shape[1]))
        self.us = np.zeros((self.grid.Nelem, self.solver.shape[1]))

    Nelem = self.grid.Nelem
    Wi = spdiags( 2*np.pi * self.grid.hr * self.grid.ha * self.grid.wi, 0, Nelem, Nelem)
    W  = spdiags( self.grid.w, 0, Nelem, Nelem )
    WWi = W @ Wi

    isolver = self.solver[:, ispin]
    Nsol = isolver.shape[0]
    Nmo, N, m = np.empty((Nsol), dtype=int), np.empty((Nsol), dtype=int), np.empty((Nsol), dtype=int)
    for i, j in enumerate(isolver):
        j.hamiltonian()
        Nmo[i] = j.Nmo
        N[i] = j.N
        m[i] = j.m
        
    #Transforming cell array from matlab to numpy array
    C = np.empty((Nsol,1), dtype=object)  
    D = np.empty((Nsol,1), dtype=object)

    for it in range(Nsol):
        # C[it] = matlab.repmat( Wi @ isolver[it].H0, Nmo[it], 1 )
        C[it,0] = np.repeat( Wi @ isolver[it].H0, Nmo[it])[:,None]

        if m[it] == 0:            
            D[it,0] = 2 * np.ones( (Nmo[it], 1) )

            if isolver[it].pol == 1:
                nu = N[it] / 2 - np.floor(N[it] / 2)
            else:
                nu = N[it] - np.floor(N[it])
            
            if nu != 0:
                D[it,0][Nmo[it]] = D[it][Nmo[it]] * nu
        
        else:
            D[it,0] = 4 *  np.ones((Nmo[it,0], 1))
            if isolver[it].pol == 1:
                nu = N[it] / 4 - np.floor(N[it] / 4)
            else:
                nu = N[it] / 2 - np.floor(N[it] / 2)

            if nu != 0:
                D[it,0][Nmo[it]] = D[it][Nmo[it]] * nu


    H_matrices = []
    for i in range(C[0,0].shape[0]):
        H_matrices.append(C[0,0][i])
    matrices = (H_matrices[i][0] for i in range(len(H_matrices)))
    H = block_diag( matrices )
    occ = np.vstack(D[:][0])

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
                         (e0[:np.sum(Nmo)-1] - e0[-1])[:, None], 
                         np.zeros((North, 1)),
                         vs0 - e0[np.sum(Nmo)-1]
                        ))

    else:
        phi   = -1.0 * isolver[0].phi
        evals = isolver[0].eig

        for i in range(1,Nsol):
            warnings.warn("WARNING, Double Check concatenation")
            phi   = np.hstack( (phi, (isolver[i].phi) ))
            evals = np.hstack( (evals, isolver[i].eig) )

        #Sanity chech. If there is one evals. It wont be able to be concatenated
        if evals.ndim == 1:
            evals = evals[:,None]

        X = np.vstack(( phi.flatten("F")[:, None],
                        (evals[0:-1] - evals[-1]),
                        np.zeros((North, 1)),
                        (isolver[0].veff - evals[-1])[:,None]
                      ))

    if self.optInv.ab_sym is True:
        X = ab_symmetrize(X, Nelem, Nmo, North, self)
    X = normalize(X, Nmo, Nelem, WWi, occ)

    if self.optInv.disp is True:
        # print('            _________________________')
        # print('           | Internal Inversion Loop |')
        # print('            _________________________')
        print('           iter      res_ks        res_ncon         res_Ncon    res_linsolve  iter_linsolve\n');
        print('           ________________________________________________________________________________\n');

    iter        = 0
    res         = 1
    finished    = False
    ForceUpdate = False

    if self.optInv.avoid_loop is True:
        old_res_ks     = 0
        old_res_ncon   = 0
        older_res_ks   = 0
        older_res_ncon = 0

    while not finished:
        iter += 1
        jac, eqn = kmatrix(X, i_orth, Nelem, Nmo, WWi, H, n0, occ, B2i, B2j, B3i, B3j)

        if self.optInv.avoid_loop is True and iter > 1:
            older_res_ks   = old_res_ks
            older_res_ncon = old_res_ncon
            old_res_ks     = res_ks
            old_res_ncon   = res_ncon

        #Inner range of phi:
        if np.sum(Nmo) - 1 == 0:
            res_Ncon = 0.0
        else:
            subset = np.array(range(np.sum(Nmo) - 1)) + np.sum(Nmo) * Nelem
            res_Ncon = np.max(np.abs( eqn[ subset ] ))

        res_ks   = np.max(np.abs( eqn[ np.array(range(Nelem * np.sum(Nmo))) ]))
        res_ncon = np.max(np.abs(
                    eqn[ np.array(range( Nelem )) + np.sum(Nmo) * Nelem + np.sum(Nmo) - 1 + North  ]
                    ))
        # res_Ncon = np.max(np.abs( eqn[ subset ] ))
        res = max(res_ks, res_ncon)

        if self.optInv.avoid_loop is True and iter > 2:
            if res      < 1e-6 and \
               res_ks   > old_res_ks     * 0.1 and \
               res_ks   < older_res_ks   * 10. and \
               res_ncon > older_res_ncon * 0.1 and \
               res_ncon < older_res_ncon * 10.:

                ForceUpdate = True
            else:
                ForceUpdate = False

        if self.optInv.disp:
             print(f"        {iter:5d}      {res_ks:7.5e}     {res_ncon:7.5e}     {res_Ncon:7.5e}")

        if res >1e2:
            print(f"      ERROR too large: {res}")
            finished = True

        if res < self.optInv.tol_invert:
            finished = True

        elif iter >= self.optInv.max_iter_invert:
            finished = True
            warnings.warn('\n Convergence not reached at maximum iteration')

        else:
            #Convergence restrictions have been met
            if (res_ncon > (res_ks * self.optInv.res_factor) or \
                res_ks   < (self.optInv.tol_invert * 1e2) or    \
               ForceUpdate) == True:

                #Use Iterative?
                if self.optInv.use_iterative == True:
                    warnings.warn("QMR solver may not work")
                    tol   = self.optInv.tol_lin_solver
                    maxit = self.optInv.max_iter_lin_solver

                    A = jac[ :(np.sum(Nmo) * Nelem), :(np.sum(Nmo) * Nelem)]
                    B = jac[ :(np.sum(Nmo) * Nelem), (np.sum(Nmo) * Nelem):]
                    M1 = np.concatenate((A,B), axis=1)
                    M2 = np.concatenate((B.T , self.optInv.k_val * eye((B.shape[1])).toarray()), axis=1)
                    M = np.concatenate((M1, M2), axis=0)
                    M = csc_matrix(M)

                    M_inv = inv(M)

                    lu = scipy.sparse.linalg.splu(M)

                    #Neet to build a Linear operator out of L and U Matrix.
                    #QMR is not written to take matrices, but Linear Operators. 
                    M1 = lu.L.A
                    M2 = lu.U.A

                    M1_solve  = lambda z: spsolve(M1, z)
                    M1_solveH = lambda z: spsolve(np.conjugate(M1).T, z)

                    M2_solve  = lambda z: spsolve(M2, z)
                    M2_solveH = lambda z: spsolve(np.conjugate(M2).T, z)

                    LO_M1 = LinearOperator(M.shape, matvec=M1_solve, rmatvec=M1_solveH)
                    LO_M2 = LinearOperator(M.shape, matvec=M2_solve, rmatvec=M2_solveH)

                    # dX, info = lgmres(jac, -eqn, M=M)
                    # # dX, info = scipy.sparse.linalg.gmres(jac, -eqn, M=M_inv)
                    # dX, info = qmr( jac, -eqn, tol=tol,maxiter=maxit, M1=LO_M1, M2=LO_M2)  
                    dX, info = qmr(jac, -eqn, M1=LO_M1, M2=LO_M2)
                    dX = dX[:,None]

                    if info > 0:
                        print(f"      Iterative solver failed | Max iterations surpased")
                        finished = True

                else:
                    dX = - 1.0 * np.linalg.solve(jac, eqn)

                #Add dX to X
                dv = dX[ np.array(range( Nelem )) + np.sum(Nmo) * Nelem + np.sum(Nmo) -1 + North ]
                X += dX

            #If conditions have not been met, recalculate X
            else:
                #Inner range of phi:
                if np.sum(Nmo) - 1 == 0:
                    subset = np.array(range(np.sum(Nmo) * Nelem))
                    evals = np.concatenate( (np.array(([[0]])) ), axis=0 )
                else:
                    subset = np.array(range(np.sum(Nmo) - 1)) + np.sum(Nmo) * Nelem
                    evals = np.concatenate( ( X[ subset ], np.array(([[0]])) ), axis=0 )

                self.vs[:, ispin] = X[np.array(range(Nelem)) + np.sum(Nmo) * Nelem + np.sum(Nmo) - 1 + North, 0]

                for it in range(Nsol):

                    isolver[it].phi = np.reshape( X[ np.array(range(Nmo[it] * Nelem)) + np.sum(Nmo[0:it-1]) * Nelem] , (Nelem, Nmo[it]) , order="F")
                    isolver[it].eig = evals[:Nmo[it] + np.sum(Nmo[:it-1])]

                for i in isolver:
                    i.setveff(self.vs[:, ispin])
                    i.calc_orbitals()

                phi = None
                evals = None
                for i in isolver:
                    phi   = i.phi if phi is None else np.concatenate( (phi, i.phi),   axis=1)
                    evals = i.eig if evals is None else np.concatenate( (evals, i.eig), axis=1)

                X = np.concatenate((phi.flatten("F")[:,None], 
                                    (evals[:-1] - evals[-1])[:, None], 
                                    np.zeros((North, 1)), 
                                    (self.vs[:, ispin]-evals[-1])[:, None]),
                                    axis=0 )
                        
            if self.optInv.ab_sym is True:
                X = ab_symmetrize(X, Nelem, Nmo, North, self)

            X = normalize(X, Nmo, Nelem, WWi, occ)

        if self.optInv.disp is True:
            pass

    #Inner range of phi:
    if np.sum(Nmo) - 1 == 0:
        subset = np.array(range(np.sum(Nmo) * Nelem))
        evals = np.concatenate( (np.array(([[0]])) ), axis=0 )
    else:
        subset = np.array(range(np.sum(Nmo) - 1)) + np.sum(Nmo) * Nelem
        evals = np.concatenate( ( X[ subset ], np.array(([[0]])) ), axis=0 )

    self.vs[:, ispin] = X[ np.array(range(Nelem)) + np.sum(Nmo) * Nelem + np.sum(Nmo) - 1 + North, 0]

    for it in range(Nsol):
        isolver[it].phi = np.reshape( X[ np.array(range(Nmo[it] * Nelem)) + np.sum(Nmo[0:it-1]) * Nelem] , (Nelem, Nmo[it]) , order="F")
        isolver[it].eig = evals[:Nmo[it] + np.sum(Nmo[:it-1])]

    for it in range(Nsol):
        self.us[ispin] = isolver[it].get_homo()

        isolver[it].calc_density()
        isolver[it].setveff(self.vs[:, ispin])

    if res > self.optInv.tol_invert:
        flag = False
    else:
        flag = True


    inv_info = inversion_info
    inv_info.nfev       = iter 
    inv_info.optimality = res
    output = np.empty((1,1), dtype=object)
    output[0,0] = inv_info

    return flag, output







