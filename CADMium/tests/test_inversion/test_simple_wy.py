import numpy as np
import CADMium as cad
import pytest

@pytest.fixture()
def he2():
    a = 4.9322/2

    Za = 2
    Zb = 2

    pol = 1

    Nmo_a = [[1]]
    Na   = [[2]]
    nua  = 1

    Nmo_b = [[1]]
    Nb    = [[2]]
    nub   = 1

    Nmo = [[2]]
    Nm  = [[4]]

    optPartition = {"AB_SYM" : False,
                    "ENS_SPIN_SYM" : True}

    NP = 7
    NM = [4, 4]
    L = np.arccosh(12/a)
    loc = np.array(range(-4,5))

    X = cad.Psgrid(NP, NM ,a, L, loc)
    X.initialize()

    ks = cad.Kohnsham(X, Za, Zb, pol, Nmo, Nm, {})
    ks.scf({})
    return ks

def test_energy(he2):
    assert(np.isclose(he2.E.E, -5.670495229536985))

def test_inversion(he2):
    #Set up Inversion
    dm = he2.n

    #Options
    optInversion = {"invert_type" : "wuyang"}

    #Set up objects
    #Partition object required for initial guess
    P  = cad.Pssolver( X, Nmo, Nm )
    WY = cad.Inverter( X, P, optInversion )
    part = cad.Partition( X, Za, Zb, pol, Nmo_a, Na, nua, Nmo_b, Nb, nub, {} )
    part.optPartition.isolated = True
    part.scf({"disp" : True, 
              "e_tol" : 1e-7})

    assert(np.isclose(part.E.E, -4.8590922158945995))

    phi0, e0, v0 = part.initialguessinvert(ispin=0) 
    success, inv_info = WY.invert( dm, v0, phi0, e0, )

    assert( success, 1 )
