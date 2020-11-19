import numpy as np
import pytest
from CADMium import Psgrid
from CADMium import Kohnsham

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
            "interaction_type" : "ni",
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


def test_kinetic(hydrogen):
    assert(np.isclose(hydrogen.E.Ts, 0.5)) 

def test_external(hydrogen):
    assert(np.isclose(hydrogen.E.Vnuc, -1.0))

def test_total(hydrogen):
    assert(np.isclose(hydrogen.E.E, -0.5))