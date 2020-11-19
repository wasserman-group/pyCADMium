import numpy as np
import pytest
from CADMium import Psgrid
from CADMium import Kohnsham

#Tested against the Atomic Reference Data for Electronic Structure Calculations from Nist
# https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-7-0

@pytest.fixture()
def lithium():
    #Distance of the nucley from grid center
    a =  1.0

    #Nuclear charges on centers AB
    Za  = 3
    Zb = 0

    #Set polaization. 1 Unpolarized, 2 Polarized
    pol = 1

    Nmo = [[2]]
    N   = [[3]]

    optKS = {
            "interaction_type" : "dft",
            "SYM" : False,
            "FRACTIONAL" : True,
            }

    #Grid Options
    NP = 7 #Number of points per block
    NM =  [4,4] #Number of blocks [angular, radial]
    L = np.arccosh(15./a) #Maximum radial coordinate value
    loc = np.array(range(-4,5)) #Non inclusive on upper bound

    #Create and initialize grid object
    grid = Psgrid(NP, NM, a, L, loc)
    grid.initialize()

    #Kohn Sham object
    KS = Kohnsham(grid, Za, Zb, pol, Nmo, N, optKS)
    KS.scf(optKS)

    return KS


def test_kinetic(lithium):
    assert(np.isclose(lithium.E.Ts, 7.235654)) 
    
def test_total(lithium):
    assert(np.isclose(lithium.E.Et, -7.335109)) 
    
def test_hartree(lithium):
    assert(np.isclose(lithium.E.Eh, 3.991242)) 
    
def test_nuclear(lithium):
    assert(np.isclose(lithium.E.Enuc, -16.909834)) 
    
def test_exc(lithium):
    assert(np.isclose(lithium.E.Ex + lithium.E.Ec, -1.491935 +  -0.160236)) 
    
def test_eigenvalue(lithium):
    assert(np.isclose(lithium.solver[0,0].eig[0], -1.87826847)) 
    assert(np.isclose(lithium.solver[0,0].eig[1], -0.105571)) 