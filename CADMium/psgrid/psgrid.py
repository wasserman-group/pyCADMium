"""

psgrid.py
Provides handling for 2d finite difference meshes on a prolate spheroidal grid

"""

import sys
sys.path.append("..")

from .initialize import initialize
from .mirror import mirror
from .square import square
from .finite_difference_1d import finite_difference_1d


class Psgrid():
    """
    Generates spheroidal grid

    Parameters
    ----------

    NP: int
        Number of points per integration block
    NM: list
        Number of angluar/radial blocks
    a: float
        half bond lenght
    L: float
        spheroidal box size
    loc:
       
    """
    def __init__(self, NP, NM, a, L, loc):

        #Mesh properties
        self.NP = NP
        self.NMa = NM[0]
        self.NMr = NM[1]
        self.Na = None
        self.Nr = None
        self.Nelem = None

        #Prolate Spheroidal Coordinates
        self.xa = None
        self.xr = None
        self.Xa = None
        self.Xr = None
        self.ha = None
        self.hr = None

        #Cartesian coordinates
        self.Y = None
        self.Z = None

        #Constants
        self.a = a
        self.R = None
        self.L = L

        #volume element and integration weights
        self.w = None
        self.wi = None
        self.f = None

        #Finite difference stencils
        self.d1 = None
        self.i1 = loc
        self.d2 = None
        self.i2 = loc

        #Basic finite difference operators
        # m -> even symmetry
        self.eDa1 = None
        self.eDa2 = None
        self.eDr1 = None
        self.eDr2 = None
        # -> odd symmetry
        self.oDa1 = None
        self.oDa2 = None
        self.oDr1 = None
        self.oDr2 = None

        #Prolate spheroidal operators
        self.elap = None
        self.olap = None
        self.grada = None
        self.gradr = None
        self.diva = None
        self.divr = None

        #Boundary Conditions
        self.bcN = None
        self.bc1 = None
        self.blap = None
        self.bXa = None
        self.bXr = None

        #Scale Factors
        self.h1 = None
        self.h2 = None
        self.h3 = None
        self.L_lap = None
        self.U_lap = None
        self.DISP = True


    #Import methods
    
    def initialize(self):
        initialize(self)

    def mirror(self, fin):
        return mirror(self, fin)


    def square(self, fin):
        return square(self, fin)

    def finite_difference_1d(self):
        finite_difference_1d(self)

    # def finite_difference_2d(self):
    #     finite_difference_2d(self)


