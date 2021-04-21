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
    if not self.ens:
        iksa = [self.KSa]
        iksb = [self.KSb]
    else:
        iksa = [self.KSa, self.KSA]
        iksb = [self.KSb, self.KSB]
    
    for KSa in iksa:
        KSa.V.v_s = Pb * np.log2( Pb/Pa ) / np.sum(self.nf, axis=1)
        KSa.V.v_s = np.repeat( KSa.V.v_s, 2, axis=1 )

    for KSb in iksb:
        KSb.V.v_s = Pa * np.log2( Pa/Pb ) / np.sum(self.nf, axis=1)
        KSb.V.v_s = np.repeat( KSb.V.v_s, 2, axis=1 )

    