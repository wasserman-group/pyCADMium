"""
scf.py
"""

import numpy as np

def scf(self, optKS):
    """
    SCF method to handle self consistent field calculations
    """

    Tolerance = optKS["Tolerance"] if "Tolerance" in optKS.keys() else 10e-6
    MaxIter = optKS["MaxIter"] if "MaxIter" in optKS.keys() else 100
    Alpha = optKS["Alpha"] if "Alpha" in optKS.keys() else 0.50
    Verbose = optKS["Verbose"] if "Verbose" in optKS.keys() else True
    CONTINUE = optKS["CONTINUE"] if "CONTINUE" in optKS.keys() else False
    ITERATIVE = optKS["ITERATIVE"] if "ITERATIVE" in optKS.keys() else True
    SPINFLIPSYM = optKS["SPINFLIPSYM"] if "SPINFLIPSYM" in optKS.keys() else False
    AutoTol = optKS["AutoTol"] if "AutoTol" in optKS.keys() else False
    AutoTolIter = optKS["AutoTolIter"] if "AutoTolIter" in optKS.keys() else 3

    if Verbose is True:
        print(' iter    Total Energy     HOMO Eigenvalue         Res       \n');
        print('----------------------------------------------------------- \n');

    if CONTINUE is True:
        #If we continue a calculation, we check that we have an input density
        assert len(self.n) == 0, "CONTINUE option is True, but there is no input density"
        self.vext = np.zeros_like(self.vnuc) if self.vext is None else self.vext

    elif CONTINUE is False:
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

    if AutoTol == True:
        min_dif = 10
        num_iter_not_min = 0


    while (diff > Tolerance or AutoTol == True and num_iter_not_min < AutoTolIter) and iter < MaxIter:


        #Calculate and set new effective potential:
        self.calc_hxc_potential() 
        self.veff = self.vnuc + self.vext + self.vhxc

        for i in range(self.Nmo.shape[0]):
            for j in range(self.Nmo.shape[1]):
                self.solver[i,j].setveff(self.veff[:,j])

        #Spin flip mirror symmetry
        if SPINFLIPSYM is True:
            self.veff = (self.veff + self.grid.mirror(self.grid.spinflip(self.veff))) / 2

        #Calculate new density
        if ITERATIVE is True and CONTINUE is True:
            nout = self.calc_density(ITERATIVE==ITERATIVE)

        else:
            nout = self.calc_density(ITERATIVE=ITERATIVE, dif=diff)


        #Set new density with linear mixing
        self.n = (1 - Alpha) * self.n + Alpha * nout

        if SPINFLIPSYM is True:
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

        if AutoTol is True:
            if dif < min_dif:
                num_iter_not_min = 0
                min_dif = diff
            else:
                num_iter_not_min = num_iter_not_min + 1

        if Verbose is True:
            print(f"   {iter}         {self.E.E:.3f}          {self.u:.3f}            {diff}")

        iter += 1


        



        