"""
paramke.py
"""
import numpy as np

import pylibxc as libxc

from ..functionals.TW_paramke import TW_paramke
from ..functionals.PW91k_paramke import PW91k_paramke

class Paramke():
    """
    Provides an interface to the libxc library
    """

    def __init__(grid, k_family, func_id, param):

        self.grid = grid
        self.k_family = k_family, 
        self.func_id = func_id
        self.param = self.param

    def get_vxc_components(self, n):

        #Calculate total energy
        pol = self.n.shape[1]
        exc = self.get_exc(n)
        Exc = self.grid.integrate(exc * np.sum(n, axis=1))

        #Calculate potential
        vxc = self.get_vxc(n)


    def get_exc(self, n):

        pol = self.n.shape[1]

        if self.family == 'gga':

            s = self.grid.reduced_grad(n)

            if self.func_id == 1001:
                #Requires as param {"kappa", "mu"}
                F = TW_paramke(s, self.param)

            elif self.func_id == 1002:
                #Requires as parameters {"A1", "A2", "A3", "A4", "A", "B1"}
                F = PW91K_paramke(s, self.param)

            else:
                raise ValueError("Funcional ID not recognized")


            C_TF = 0.3 * (3 * np.pi ** 2) ** (2/3)

            if pol == 1:
                exc = 0.5 * C_TF @ n**(2/3) * F
                #Missing warning to avoid infty in exc#
            elif pol == 2:
                zk = 0.5 * C_TF @ (2 * n)**(5/3) * F
                exc = np.sum(zk,axis=1) / np.sum(n, axis=1)
        
        else: 
            raise ValueError("Only gga family available")

    def get_vxc(self, n):

        pol = n.shape[1]
        #Potential contribution

        if self.k_family == "lda":
            raise ValueError("There is no support for LDA family paramke.py")

        elif self.k_family == "gga":
            sigma = self.grid.sigma(n)

            vrho, vsigma = gga_vxc(n, sigma, pol)

            raise ValueError("Not Yet Implemented")
        
        
    def gga_vxc(self, n, sigma, pol):

        sval = {"f" : s, "ds": np.ones((n.shape[0], n.shape[1]))}

        if self.func_id == 1001:
            Raise ValueError("Automatic Differenciation not yet Implemented")
            #Required parameters are {"kappa", "mu"}
            #F = TW_paramke(s, self.param)
            #dFds = dTW_paramke(sval, self.param)

        if self.func_id == 1002
            Raise ValueError("Automatic Differenciation not yet implemented")


        C_TF  = 0.3 * (3*np.pi**2)**(2/3)

        if pol == 1:
            dsdn = -2 * @ sigma ** (0.5) / 3 / (3*np.pi**2)**(1/3) / n**(7/3)
            dsdsimga = 1/4/(3*np.pi**2)**(1/3) / sigma**0.5 / n**(4/3)

            #DFDS may not be given right now
            vsigma =  C_TF @ n**(5/3) * dFds["ds"] * dsdsigma
            vrho   = 5/3 * C_TF * n**(2/3) * F + C_TF * n**(5/3) * dFds["ds"] * dsdn

        if pol == 2:
            raise ValueError("Not Implemented Yet")