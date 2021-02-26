"""
vp_overalp.py
Do calculations for overlap type functionals
"""

import numpy as np
from scipy.special import erf

def vp_overlap(self):

    const = 2

    #Calculate overlap
    self.E.S = self.grid.integrate((np.sum(self.na_frac, axis=1) * np.sum(self.nb_frac, axis=1))**(0.5))
    self.E.F = erf( const * self.E.S)

    #Functional derivative of the overlap
    self.KSa.V.dSdn = self.KSa.scale * 0.5 * (self.nb_frac / self.na_frac)**0.5
    if self.optPartition is True:
        self.KSa.V.dSdn = np.repeat(self.KSa.V.dSdn, 2, axis=1)
    #Remove any Nans
    self.KSa.V.dSdn[ np.isinf(self.KSa.V.dSdn) ] = 0.0
    self.KSa.V.dSdn[ np.isnan(self.KSa.V.dSdn) ] = 0.0
    self.KSa.V.dFdn = 2 * np.pi**(-0.5) * np.exp( -(const * self.E.S)**2 ) * const * self.KSa.V.dSdn 

    #Functional derivative of the overlap
    self.KSb.V.dSdn = self.KSb.scale * 0.5 *  (self.na_frac / self.nb_frac)**0.5
    if self.optPartition is True:
        self.KSb.V.dSdn = np.repeat(self.KSb.V.dSdn, 2, axis=1)
    #Remove any Nans
    self.KSb.V.dSdn[ np.isinf(self.KSb.V.dSdn) ] = 0.0
    self.KSb.V.dSdn[ np.isnan(self.KSb.V.dSdn) ] = 0.0
    self.KSb.V.dFdn = 2 * np.pi**(-0.5) * np.exp( -(const * self.E.S)**2 ) * const * self.KSb.V.dSdn 


