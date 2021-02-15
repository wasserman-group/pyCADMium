import numpy as np
import pytest
from CADMium import Psgrid
from CADMium import Kohnsham

#Tested against the Atomic Reference Data for Electronic Structure Calculations from Nist
# https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-7-0

@pytest.fixture()
def hydrogen():
    #Distance of the nucley from grid center
    a =  1.0

    #Nuclear charges on centers AB
    Za  = 1
    Zb = 0

    #Set polaization. 1 Unpolarized, 2 Polarized
    pol = 1

    Nmo = [[1]]
    N   = [[1]]

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
    KS.scf()

    return KS


def test_kinetic(hydrogen):
    print(hydrogen.E.Ts)
    assert(np.isclose(hydrogen.E.Ts, 0.424864)) 
    
def test_total(hydrogen):
    assert(np.isclose(hydrogen.E.Et, -0.445667)) 
    
def test_hartree(hydrogen):
    assert(np.isclose(hydrogen.E.Eh, 0.282770)) 
    
def test_nuclear(hydrogen):
    assert(np.isclose(hydrogen.E.Enuc, -0.920823)) 
    
def test_exc(hydrogen):
    assert(np.isclose(hydrogen.E.Ex + hydrogen.E.Ec, -0.193072 + -0.039407)) 
    
def test_eigenvalue(hydrogen):
    assert(np.isclose(hydrogen.solver[0,0].eig[0], -0.233454)) 