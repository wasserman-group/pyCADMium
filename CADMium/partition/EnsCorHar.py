"""
EnsCorHar.py
Calculates ensemble correlated hartree energy
"""

import numpy as np

def EnsCorHar(self):
    if True:
        print("Warning: If len(KS) > 1 Has not been migrated from matlab")
        if self.optPartition.ens_spin_sym is not True:
            self.KSa.E.Ehcor = 0.0
            self.KSb.E.Ehcor = 0.0
            self.KSa.V.vhcor = np.zeros_like( self.KSa.V.vh )
            self.KSb.V.vhcor = np.zeros_like( self.KSb.V.vh )
            self.E.Ehcor = 0.0

        else:
            self.KSa.E.Ehcor = np.sum( ( self.grid.integrate(self.KSa.n * self.grid.spinflip(self.KSb.V.vh) ) ) ) / 2
            self.KSb.E.Ehcor = np.sum( ( self.grid.integrate(self.KSb.n * self.grid.spinflip(self.KSa.V.vh) ) ) ) / 2
            
            self.KSa.V.vhcor = self.grid.spinflip( self.KSb.V.vh )
            self.KSb.V.vhcor = self.grid.spinflip( self.KSa.V.vh )

            self.E.Ehcor = 0.0
            for iks in [self.KSa, self.KSb]:
                self.E.Ecor = self.E.Ehcor + iks.scale * 2 * iks.E.Ehcor


        
