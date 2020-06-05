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
from .finite_difference_2d import finite_difference_2d
from .operators import operators


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
        self.NP = NP        #Number of points per integration block
        self.NMa = NM[0]    #Number of angular blocks
        self.NMr = NM[1]    #Number of radial blocks
        self.Na = None      #Number of angular points
        self.Nr = None      #Number of radial points
        self.Nelem = None   #Total number of elements (points?)

        #Prolate Spheroidal Coordinates
        self.xa = None      #Angular coordinate in 1d array
        self.xr = None      #Radial coordinate in 1d array
        self.Xa = None      #Angular coordinate in 2d grid
        self.Xr = None      #Radial coordinate in 2d grid
        self.ha = None      #Angular grid spacing
        self.hr = None      #Radial grid spacing

        #Cartesian coordinates
        self.Y = None
        self.Z = None

        #Constants
        self.a = a          #Half bond length
        self.R = None       #Bond length
        self.L = L          #Spheroidal box size

        #volume element and integration weights
        self.w = None       #Volume element
        self.wi = None      #Integration weights
        self.f = None       #Orbital angular momentum potential

        #Finite difference stencils
        self.d1 = None      #First order coefficients
        self.i1 = loc       #Location
        self.d2 = None      #Second order coefficients
        self.i2 = loc       #Location

        #Basic finite difference operators
        # m -> even symmetry
        self.eDa1 = None    #Angular differentiator
        self.eDa2 = None    #Angular differentiator
        self.eDr1 = None    #Radial differentiator
        self.eDr2 = None    #Radial differentiator
        # -> odd symmetry
        self.oDa1 = None    #Angular differentiator
        self.oDa2 = None    #Angular differentiator
        self.oDr1 = None    #Radial differentiator
        self.oDr2 = None    #Radial differentiator

        #Prolate spheroidal operators
        self.elap = None    #Laplacian -> even
        self.olap = None    #Laplacian -> odd
        self.grada = None   #Angular gradient component (we only need m=even gradient)
        self.gradr = None   #Radial gradient component
        self.diva = None    #Angular divergence component (we only need m=even gradient)
        self.divr = None    #Radial divergence component

        #Boundary Conditions
        self.bcN = None     #Size of boundary region
        self.bc1 = None     #Outer radial boundary conditions 1st order
        self.bc2 = None     #Outer radial boundary conditions 2nd order
        self.blap = None    #Laplacian for values beyond Xr=L boundary
        self.bXa = None     #Coordinates just outsise the Xr=Boundary
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

    def finite_difference_2d(self):
        finite_difference_2d(self)

    def operators(self):
        operators(self)


