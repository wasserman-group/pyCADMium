"""
scf.py
"""

import numpy as np
from pydantic import BaseModel

class KohnShamSCFOptions(BaseModel):
    disp          : bool = True
    continuing    : bool = False
    iterative     : bool = True
    spinflipsym   : bool = False
    auto_tol       : bool = False
    maxiter       : int = 100
    auto_tol_iter : int = 3
    e_tol         : float = 10e-6
    alpha         : float = 0.5



def scf(self, optKS):
    """
    SCF method to handle self consistent field calculations
    """

    for i in optKS.keys():
        i = i.lower()
        if i not in KohnShamSCFOptions().dict().keys():
            raise ValueError(f"{i} is not a valid option for the SCF procedure")
    optKS = KohnShamSCFOptions(**optKS)

    if optKS.disp is True:
        print(' iter    Total Energy     HOMO Eigenvalue         Res       \n');
        print('----------------------------------------------------------- \n');

    if optKS.continuing is True:
        print("Am I continuing?")
        #If we continue a calculation, we check that we have an input density
        assert len(self.n) != 0, "CONTINUE option is True, but there is no input density"
        self.vext = np.zeros_like(self.vnuc) if self.vext is None else self.vext
    else:
        #We need an initial guess
        self.vext = np.zeros_like(self.vnuc) if self.vext is None else self.vext
        self.vhxc = np.zeros_like(self.vnuc)
        #Initial guess for effective potential is just nuclear potential
        self.set_effective_potential()

        nout = self.calc_density()

        #Set the new input density
        self.n = nout

    #Start up the scf loop
    diff  = 10
    old_E = 0
    old_n = np.zeros((self.grid.Nelem, self.pol))
    iter = 1

    if optKS.auto_tol == True:
        min_dif = 10
        num_iter_not_min = 0


    while (diff > optKS.e_tol or optKS.auto_tol == True and num_iter_not_min < opt.KS.auto_tol_iter) and iter < optKS.maxiter:


        #Calculate and set new effective potential:
        self.calc_hxc_potential() 
        self.veff = self.vnuc + self.vext + self.vhxc

        for i in range(self.Nmo.shape[0]):
            for j in range(self.Nmo.shape[1]):
                self.solver[i,j].setveff(self.veff[:,j])

        #Spin flip mirror symmetry
        if optKS.spinflipsym is True:
            self.veff = (self.veff + self.grid.mirror(self.grid.spinflip(self.veff))) / 2

        #Calculate new density
        if optKS.iterative is True and optKS.continuing is True:
            nout = self.calc_density(ITERATIVE=optKS.iterative)
        else:
            nout = self.calc_density(ITERATIVE=optKS.iterative, dif=diff)


        #Set new density with linear mixing
        self.n = (1 - optKS.alpha) * self.n + optKS.alpha * nout

        if optKS.spinflipsym is True:
            self.n = (self.n + self.grid.mirror(self.grid.spinflip(self.n)))/2
            

        #Calculate energies and chemical potential
        self.energy()
        self.calc_chempot()

        #Convergence check
        dif_E = np.abs( (self.E.E - old_E) / self.E.E )
        dif_n = np.max(  self.grid.integrate(np.abs(self.n - old_n)) / self.grid.integrate(np.abs(self.n))  )


        diff = max(dif_E, dif_n)
        old_E = self.E.E
        old_n = self.n

        if optKS.auto_tol is True:
            if dif < min_dif:
                num_iter_not_min = 0
                min_dif = diff
            else:
                num_iter_not_min = num_iter_not_min + 1

        if optKS.disp is True:
            print(f"{iter:5d}      {self.E.E:+9.5f}      {self.u:+9.5e}       {diff:+9.5e}")

        iter += 1
