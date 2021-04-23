"""
EnsCorHar.py
Calculates ensemble correlated hartree energy
"""

import numpy as np

def EnsCorHar(self):
    if not self.ens:  # No Ensemble
        if not self.optPartiton.ens_spin_sym: #Not ENS Spin Symmetry
            self.KSa.E.Ehcor = 0.0
            self.KSb.E.Ehcor = 0.0
            self.KSa.V.vhcor = np.zeros_like( self.KSa.V.vh )
            self.KSb.V.vhcor = np.zeros_like( self.KSb.V.vh )
            self.E.Ehcor = 0.0
        else:                                 #With ENS Spin Symmetry
            self.KSa.E.Ehcor = self.grid.integrate( np.sum(self.KSa.n * self.grid.spinflip(self.KSb.V.vh) ))   / 2.0
            self.KSb.E.Ehcor = self.grid.integrate( np.sum(self.KSb.n * self.grid.spinflip(self.KSa.V.vh) ))   / 2.0
            self.KSa.V.vhcor = self.grid.spinflip( self.KSb.V.vh )
            self.KSb.V.vhcor = self.grid.spinflip( self.KSa.V.vh )
            self.E.Ehcor = 0.0
            for iks in [self.KSa, self.KSb]:
                self.E.Ehcor += iks.scale * 2 * iks.E.Ehcor
    else:             # Yes Ensemble

        # print("Potencial de Hartree\n")

        # print("A ver cuanto sale de esto?", self.grid.integrate( np.sum(self.KSa.n * self.KSB.V.vh, axis=1))/ 2  ) 



        self.KSa.E.Ehcor = self.grid.integrate( np.sum(self.KSa.n * self.KSB.V.vh, axis=1) ) / 2.0
        self.KSA.E.Ehcor = self.grid.integrate( np.sum(self.KSA.n * self.KSb.V.vh, axis=1) ) / 2.0
        self.KSb.E.Ehcor = self.grid.integrate( np.sum(self.KSb.n * self.KSA.V.vh, axis=1) ) / 2.0
        self.KSB.E.Ehcor = self.grid.integrate( np.sum(self.KSB.n * self.KSa.V.vh, axis=1) ) / 2.0

        self.KSa.V.vhcor = self.KSB.V.vh
        self.KSA.V.vhcor = self.KSb.V.vh
        self.KSb.V.vhcor = self.KSA.V.vh
        self.KSB.V.vhcor = self.KSa.V.vh

        self.E.Ehcor = 0.0
        for iks in [self.KSa, self.KSA, self.KSb, self.KSB]:
            self.E.Ehcor += iks.scale * iks.E.Ehcor       
