"""

psgrid.py
Provides handling for 2d finite difference meshes on a prolate spheroidal grid

"""

import sys
sys.path.append("..")

import numpy as np

from .initialize import initialize
from .mirror import mirror
from .spinflip import spinflip
from .square import square
from .sigma import sigma
from .integrate import integrate
from .finite_difference_1d import finite_difference_1d
from .finite_difference_2d import finite_difference_2d
from .operators import operators
from .factorize_laplacian import factorize_laplacian
from .reduced_grad import reduced_grad
from .plotter import plotter

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
    loc: np.ndarray
        stencil required for derivatives
       
    Attributes
    ----------
    Na: int
        Number of angular points
    Nr: int
        Number or radial points
    Nelem: int
        Total number of points
    xa : np.ndarray
        Angular coordinate
    xr : np.ndarray
        Radial coordinate
    Xa : np.ndarray
        Angular coordinate in 2D grid
    Xr : np.ndarray
        Angular coordiante in 2D grid
    ha : float
        Angular grid spacing
    hr : float
        Radial grid spacing    
    Y : np.ndarray
        Y axis cartesian representation of PS grid
    Z : np.ndarray
        Z axis cartesian representation of PS grid
    a : float
        Half bond length
    R : float
        Bond length
    L : float
        Spheroidal box size
    w : np.ndarray
        Volume element
    wi : np.ndarray
        Integration weights
    f : np.ndarray
        Orbital angular momentum poential
    d1 : np.ndarray
        First order coefficients
    i1 : np.ndarray
        Location of coefficients
    d2 : np.ndarray
        Second order coefficients
    i2 : np.ndarray
        Location of coefficients
    eDa1 : np.ndarray
        Angular Differentiator (Even symmetry)
    eDa2 : np.ndarray
        Angular Differentiator (Even symmetry)
    eDr1 : np.ndarray
        Radial Differentiator (Even symmetry)
    eDr2 : np.ndarray
        Radial Differentiator (Even symmetry)
    oDa1 : np.ndarray
        Angular Differentiator (Odd symmetry)
    oDa2 : np.ndarray
        Angular Differentiator (Odd symmetry)
    oDr1 : np.ndarray
        Radial Differentiator (Odd symmetry)
    oDr2 : np.ndarray
        Radial Differentiator (Odd symmetry)
    elap : csc_matrix
        Laplacian -> Even
    olap : csc_matrix
        Laplacian -> Odd
    grada : csc_matrix
        Angular gradient component
        (We only need m+=even gradient)
    gradr : csc_matrix
        Radial gradient component
    diva :  csc_matrix
        Angular divergence component
        (We only need m=even gradient)
    divr : csc_matrix
        Radial divergence component
    bcN : int
        Size of boundary region
    bc1 : np.ndarray
        Outer radial boundary conditions 1st order
    bc2 : np.ndarray
        Outer radial boundary conditions 2nd order
    blap : csc_matrix
        Laplacian for balues beond Xr=L boundary
    bXa : np.ndarray
        Coordinates just outside the Xr=L boundary
    bXr : np.ndarray
    h1 : np.ndarray
    h2 : np.ndarray
    h3 : np.ndarray
    L_lap : csc_matrix
    U_lap : csc_matrix
    DISP : logical
        Displays information about current run

    Methods
    ----------
    initialize()
        Initalizes prolate spheroidal grid
    mirror(fin)
        Mirror function accros AB plane
    square(fin)
    sigma(fin)
        Calculates gradient squared
    spinflip(fin)
        Flip spins
    integrate(f)
        Integrates a function f
    finite_difference_1d()
        Build finite difference operator matrices
    finite_difference_2d()
        Build finite difference operator matrices
    operators()
        Construct PS operators
    factorize_laplacian(DISP)
        Factorizes Laplacian for Hartree calculation
    reduced_grad(n)
        Calculates the reduced density gradient
    plotter(fin, max=1, sym=1)
        Plots function of psgrid
    """

    def __init__(self, NP, NM, a, L, loc):

        #Mesh properties
        self.NP = NP        #Number of points per integration block
        self.NMa = NM[0]    #Number of angular blocks
        self.NMr = NM[1]    #Number of radial blocks
        self.Na = None      #Number of angular points
        self.Nr = None      #Number of radial points
        self.Nelem = None   #Total number of points

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

    def sigma(self, n):
        return sigma(self, n)

    def spinflip(self, fin):
        return spinflip(self, fin)

    def integrate(self, f):
        return integrate(self, f)

    def finite_difference_1d(self):
        finite_difference_1d(self)

    def finite_difference_2d(self):
        finite_difference_2d(self)

    def operators(self):
        operators(self)

    def factorize_laplacian(self, DISP):
        factorize_laplacian(self, DISP)

    def reduced_grad(self, n):
        reduced_grad(self, n)

    def plotter(self, fin, max=1, sym=1):
        full, z, x = plotter(self, fin, max, sym)
        return full, z, x

    def plot_along_axis(self, fin, max =1, sym=1):
        f, x, y = plotter(self, fin, max, sym)
        if np.mod(f.shape[1], 2) == 1:
            mid = int(f.shape[0]/2) + 1
        else: 
            mid = int(f.shape[0]/2)
        
        return y[mid,:], f[mid,:]



