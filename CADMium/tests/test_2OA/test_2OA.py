import numpy as np
import matplotlib.pyplot as plt
from CADMium import Pssolver, Psgrid, Partition, Inverter
import pytest


@pytest.fixture
def part():
    a = 2
    #Nuclear charge for fragments A and B
    Za, Zb = 2,2
    #Set polarization 1-Unpolarized, 2-Polarized
    pol = 1
    #Fragment a electrons [alpha, beta]
    Nmo_a = [[1]] #Number of molecular orbitals to calculate
    N_a   = [[2]]
    #Ensemble mix
    nu_a = 1
    #Fragment b electrons
    Nmo_b = [[1]]
    N_b   = [[2]]
    #Ensemble mix
    nu_b = 1

    #Molecular elctron configuration
    Nmo_m = [[2]]
    N_m   = [[4]]


    #Set up grid
    NP = 2
    NM = [32,32]
    L = np.arccosh(14/a)
    loc = np.array(range(-4,5)) #Stencil outline

    grid = Psgrid(NP, NM, a, L, loc)
    grid.initialize()

    part = Partition(grid, Za, Zb, pol, Nmo_a, N_a, nu_a, Nmo_b, N_b, nu_b, { "kinetic_part_type" : "twoorbital",
                                                                            "ab_sym"            : False,
                                                                            "ens_spin_sym"      : False})
    #Setup inverter object
    mol_solver = Pssolver(grid, Nmo_m, N_m, {"tol_orbital" : 1e-9})
    part.inverter = Inverter(grid, mol_solver, {"invert_type"    : "wuyang", 
                                                "ab_sym"         : False,
                                                "ens_spin_sym"   : False,
                                                "tol_lin_solver" : 1e-3,
                                                "tol_invert"     : 1e-4,
                                                "res_factor"     : 0,
                                            })

    part.optPartition.isolated   = False
    part.scf({"disp"       : True,
            "alpha"      : [0.3],
            "max_iter"   : 1,
            "e_tol"      : 1e-7,
            "continuing" : False, 
            "iterative"  : False})

    return part

def test_vt(part):

    dT_dn1 = np.array( [-3.25473744e+02, -2.33087008e+01,  1.71899213e-01,  3.29167330e+00,
                        -3.23054983e+02, -2.53695038e+01, -1.28303564e-01,  2.94625323e-01,
                        -3.22724309e+02, -2.65481060e+01, -3.47417374e-01, -2.53168552e-01,
                        -3.20063876e+02, -2.44140930e+01,  2.09617763e+00,  2.20983076e+00])

    dT_dn2 = np.array( [ 2.41298913, -0.68025993,  2.65534717, -1.89118932, -0.58405885,
                        -0.98046271,  0.59454415,  0.5275719 , -1.13185273, -1.19957652,
                        -0.58405804,  0.85824573,  1.33114658,  1.24401849,  1.54995501,
                        3.5186795 ] )


    assert np.isclose( part.V.vt[:,0].all(), dT_dn1.all() )
    assert np.isclose( part.V.vt[:,1].all(), dT_dn2.all() )