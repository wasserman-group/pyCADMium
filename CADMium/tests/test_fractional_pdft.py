import numpy as np
import pytest
from CADMium import Pssolver
from CADMium import Psgrid
from CADMium import Partition
from CADMium import Inverter

@pytest.fixture()
def h2plus():
    a = 2.0/2
    #Nuclear charge for fragments A and B
    Za, Zb = 1,1
    #Set polarization 1-Unpolarized, 2-Polarized|
    pol = 2
    #Fragment a electrons [alpha, beta]
    Nmo_a = [[1  ,0]] #Number of molecular orbitals to calculate
    N_a   = [[0.5,0]]
    #Ensemble mix
    nu_a = 1
    #Fragment b electrons
    Nmo_b = [[1  ,0]]
    N_b   = [[0.5,0]]
    #Ensemble mix
    nu_b = 1

    #Molecular elctron configuration
    Nmo_m = [[1,1]]
    N_m   = [[1,1]]

    #Set up grid
    NP = 7
    NM = [4,4]
    L = np.arccosh(10/a)
    loc = np.array(range(-4,5)) #Stencil outline

    grid = Psgrid(NP, NM, a, L, loc)
    grid.initialize()

    part = Partition(grid, Za, Zb, pol, Nmo_a, N_a, nu_a, Nmo_b, N_b, nu_b, {  "ab_sym"            : True,
                                                                            "ens_spin_sym"      : False,  
                                                                            "kinetic_part_type" : "libxcke",
                                                                            "k_family"          : "gga",
                                                                            "ke_func_id"        : 500,
                                                                            "interaction_type"  : "ni",
                                                                            "fractional"        : True,
                                                                                })

    #Setup inverter object
    mol_solver = Pssolver(grid, Nmo_m, N_m)
    part.inverter = Inverter(grid, mol_solver, {  "AB_SYM"         : True,
                                                "ENS_SPIN_SYM"   : False,  
                                                "use_iterative"  : False,
                                                "invert_type"    : "orbitalinvert",
                                                "Tol_lin_solver" : 1e-3,
                                                "disp"           : True,  
                                                })

    # Isolated Fragment Calculation
    part.optPartition.isolated = True
    part.scf({"disp"     : False,
              "e_tol"    : 1e-7})

    part.optPartition.isolated = False
    part.scf({'iterative' : False,
                'disp'      : False,
                'continuing' : True})

    return part


def test_frag_energy(h2plus):
    assert(np.isclose(h2plus.E.Ea, -0.22348683901197508)) 

def test_mol_energy(h2plus):
    assert(np.isclose(h2plus.E.Ef, -0.44697367802395016)) 

def test_kinetic_fragment_energy(h2plus):
    assert(np.isclose(h2plus.E.Tsf, 0.6872831483178004)) 

def test_frag_KS_energy(h2plus):
    assert(np.isclose(h2plus.E.Eksf[0][0], -1.10263421)) 

def test_frag_nuclear(h2plus):
    assert(np.isclose(h2plus.E.Enucf, -1.1342568263417505)) 

def test_partition_energy(h2plus):
    assert(np.isclose(h2plus.E.Ep, -0.6556597136039474)) 

def test_ep_energy(h2plus):
    assert(np.isclose(h2plus.E.Ep_pot, -0.5704371578257259)) 

def test_ep_kinetic_energy(h2plus):
    assert(np.isclose(h2plus.E.Ep_kin, -0.08522255577822146)) 

def test_e_hxc(h2plus):
    assert(np.isclose(h2plus.E.Ep_hxc, 0)) 

def test_e_kinetic(h2plus):
    assert(np.isclose(h2plus.E.Et, -1.1026333916278976)) 

def test_total_energy(h2plus):
    assert(np.isclose(h2plus.E.E, -0.6026333916278976)) 




