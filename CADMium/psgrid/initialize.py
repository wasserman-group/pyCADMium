"""

initialize.py
Initializes grid

"""
import sys
sys.path.append("..")

import numpy as np

from ..common.NC import NC

def initialize(self, DISP=False, factor=True):

    #Number of points 
    self.Na = (self.NP-1) * self.NMa + 1 
    self.Nr = (self.NP-1) * self.NMr + 1
    self.Nelem = self.Na * self.Nr 

    #Coordinate properties
    self.xa = np.linspace(0, np.pi, self.Na+1)
    self.ha = self.xa[1] - self.xa[0]
    self.xa = self.xa + self.ha/2
    self.xa = self.xa[:-1]
    self.xa = (self.xa - np.flipud(self.xa - np.pi))/2     

    self.xr = np.linspace(0, self.L, self.Nr+1)
    self.hr = self.xr[1] - self.xr[0]
    self.xr = self.xr + self.hr/2
    self.xr = self.xr[:-1]

    Nrange = len(self.i1)

    #Boundary region
    self.bcN = int(np.floor(max([self.i1.shape[0], self.i2.shape[0]])/2))
    self.bXr = np.linspace(1, self.bcN,self.bcN)*self.hr + self.xr[-1]
    self.bXr, self.bXa = np.meshgrid(self.bXr, self.xa) 
    self.bXa = np.reshape(self.bXa, self.Na * self.bcN, order="F")
    self.bXr = np.reshape(self.bXr, self.Na * self.bcN, order="F")

    #Second grid
    self.Xr, self.Xa = np.meshgrid(self.xr, self.xa)
    self.Xa = np.reshape(self.Xa, self.Na * self.Nr, order="F")
    self.Xr = np.reshape(self.Xr, self.Na * self.Nr, order="F")

    #Cartesian coordinates
    self.Z = self.a * np.cosh(self.Xr) * np.cos(self.Xa)
    self.Y = self.a * np.sinh(self.Xr) * np.sin(self.Xa)
    self.Z = 1./2. * (self.Z - self.mirror(self.Z))
    self.Y = 1./2. * (self.Y + self.mirror(self.Y))

    if DISP is True:
        print(" Constructing integration weights ... \n")

    #Construct integration weights 1d grid
    wa = NC(self.NP, self.NMa) 
    wr = NC(self.NP, self.NMr) 
    wa = wa.reshape(len(wa), 1)
    wr = wr.reshape(len(wr), 1)
    
    #2d Grid of integration weights
    self.wi = np.reshape(wa@wr.T, (self.Nelem))
    
    if DISP is True:
        print(" Building finite difference opperators ... \n")

    #Build finite difference operator matrices
    self.finite_difference_1d()
    self.finite_difference_2d()

    #Construct prolate spheroidal operators
    self.operators()

    #Volume element
    self.w = self.a**3 * np.sin(self.Xa) * np.sinh(self.Xr) * (np.sin(self.Xa)**2 + np.sinh(self.Xr)**2)
    self.w = 0.5 * (self.w + self.mirror(self.w))

    #Angular momentum potential
    self.f = (2. / (2*self.a)**2 ) * 1./(np.cosh(self.Xr)**2 - np.cos(self.Xa)**2) \
             * (1. / (np.cosh(self.Xr)**2 - 1) + 1./(1-np.cos(self.Xa)**2))
    self.f = 0.5 * (self.f + self.mirror(self.f))

    #Factorize 
    if factor is True:
        self.factorize_laplacian(DISP)





