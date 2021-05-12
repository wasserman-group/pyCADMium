"""
vp_surpise.py
Do calculations for surprisal type functionals
"""

import numpy as np

def vp_surprise(self):

    #Calculate overlap
    Pa = np.sum(self.na_frac, axis=1) / np.sum(self.nf, axis=1)
    Pb = np.sum(self.nb_frac, axis=1) / np.sum(self.nf, axis=1)

    self.V.s = -Pa * np.log2(Pa) - Pb * np.log2(Pb)

    #Functional Derivative
    self.KSa.V.v_s = Pb * np.log2( Pb/Pa ) / np.sum(self.nf, axis=1)
    self.KSa.V.v_s = np.repeat( self.KSa.V.v_s, 2, axis=1 )

    self.KSb.V.v_s = Pa * np.log2( Pa/Pb ) / np.sum(self.nf, axis=1)
    self.KSb.V.v_s = np.repeat( self.KSb.V.v_s, 2, axis=1 )

    