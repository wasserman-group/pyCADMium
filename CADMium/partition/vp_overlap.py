"""
vp_overlap.py
Do calculations for overlap type functionals
"""

import numpy as np
from scipy.special import erf

def vp_overlap(self):

    const = 2

    #Calculate overlap
    self.E.S = self.grid.integrate((np.sum(self.na_frac, axis=1) * np.sum(self.nb_frac, axis=1))**(0.5))
    self.E.F = erf( const * self.E.S)

    if not self.ens:
        iksa = [self.KSa]
        iksb = [self.KSb]
    else:
        iksa = [self.KSa, self.KSA]
        iksb = [self.KSb, self.KSB]
    
    for KSa in iksa:
        #Functional derivative of the overlap
        KSa.V.dSdn = KSa.scale * 0.5 * (self.nb_frac / self.na_frac)**0.5
        if self.optPartition is True:
            KSa.V.dSdn = np.repeat(KSa.V.dSdn, 2, axis=1)
        #Remove any Nans
        KSa.V.dSdn[ np.isinf(KSa.V.dSdn) ] = 0.0
        KSa.V.dSdn[ np.isnan(KSa.V.dSdn) ] = 0.0
        KSa.V.dFdn = 2 * np.pi**(-0.5) * np.exp( -(const * self.E.S)**2 ) * const * KSa.V.dSdn 

    for KSb in iksb:
        #Functional derivative of the overlap
        KSb.V.dSdn = KSb.scale * 0.5 *  (self.na_frac / self.nb_frac)**0.5
        if self.optPartition is True:
            KSb.V.dSdn = np.repeat(KSb.V.dSdn, 2, axis=1)
        #Remove any Nans
        KSb.V.dSdn[ np.isinf(KSb.V.dSdn) ] = 0.0
        KSb.V.dSdn[ np.isnan(KSb.V.dSdn) ] = 0.0
        KSb.V.dFdn = 2 * np.pi**(-0.5) * np.exp( -(const * self.E.S)**2 ) * const * KSb.V.dSdn 


