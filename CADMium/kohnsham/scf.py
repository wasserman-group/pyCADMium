"""
scf.py
"""

import numpy as np

def scf(self, optKS):
    """
    SCF method to handle self consistent field calculations
    """

    Tolerance = optKS["Tolerance"] if "Tolerance" in optKS.keys() else 10e-7
    MaxIter = optKS["MaxIter"] if "MaxIter" in optKS.keys() else 50
    Alpha = optKS["Alpha"] if "Alpha" in optKS.keys() else 0.82
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

        