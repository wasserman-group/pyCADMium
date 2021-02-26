import numpy as np
from CADMium import Pssolver, Psgrid, Partition, Inverter


def test_partition_energies():

    a = 4.522/2
    #Nuclear charge for fragments A and B
    Za, Zb = 4,4
    #Set polarization 1-Unpolarized, 2-Polarized
    pol = 1
    #Fragment a electrons [alpha, beta]
    Nmo_a = [[2]] #Number of molecular orbitals to calculate
    N_a   = [[4]]
    #Ensemble mix
    nu_a = 1
    #Fragment b electrons
    Nmo_b = [[2]]
    N_b   = [[4]]
    #Ensemble mix
    nu_b = 1

    #Molecular elctron configuration
    Nmo_m = [[4]]
    N_m   = [[8]]


    #Set up grid
    NP = 4
    NM = [3,3]
    L = np.arccosh(15/a)
    loc = np.array(range(-4,5)) #Stencil outline

    grid = Psgrid(NP, NM, a, L, loc)
    grid.initialize()

    part = Partition(grid, Za, Zb, pol, Nmo_a, N_a, nu_a, Nmo_b, N_b, nu_b, { "kinetic_part_type" : "inversion",
                                                                            "ab_sym"            : True,
                                                                            "ens_spin_sym"      : False})
    #Setup inverter object
    mol_solver = Pssolver(grid, Nmo_m, N_m, {"tol_orbital" : 1e-9})
    part.inverter = Inverter(grid, mol_solver, {"invert_type"    : "wuyang", 
                                                "disp"           : True,
                                                "ab_sym"         : True,
                                                "ens_spin_sym"   : False,
                                                "tol_lin_solver" : 1e-3,
                                                "tol_invert"     : 1e-4,
                                                "res_factor"     : 0,
                                            })


    part.optPartition.isolated = True
    part.scf({"disp"  : False,
            "alpha" : [0.6],
            "e_tol" : 1e-7})

    part.optPartition.isolated   = False
    part.scf({"disp"       : True,
            "alpha"      : [0.3],
            "max_iter"   : 200,
            "e_tol"      : 1e-7,
            "continuing" : True, 
            "iterative"  : False})

    expected = {'Ea': -14.620725296666407,
                'Eb': -14.620725296666407,
                'Ef': -29.241450593332814,
                'Tsf': 23.94674969252165,
                'Eksf': np.array([[-16.43163514]]),
                'Enucf': -62.51175647729054,
                'Exf': -4.41266533612583,
                'Ecf': -0.44288284078384776,
                'Ehf': 14.179104368345756,
                'Vhxcf': 21.974617745569216,
                'Ep': -3.538585137136643,
                'Ep_pot': -7.103238227474018,
                'Ep_kin': 0.06491019980946788,
                'Ep_hxc': 3.499742890527907,
                'Et': -32.78003573046946,
                'Vnn': 3.5382574082264484,
                'E': -29.241778322243007,
                'evals_a': np.array([], dtype=np.float64),
                'evals_b': np.array([], dtype=np.float64),
                'Ep_h': 3.548294203236745,
                'Ep_x': -0.04334766216506569,
                'Ep_c': -0.0052036505437722536}

    for i in part.E.__dict__:
        if i.startswith("__") is False:
            if type(getattr(part.E, i)) == np.ndarray:
                assert np.isclose(  getattr(part.E, i).all(), expected[i].all()  )
            else:
                assert np.isclose(  getattr(part.E, i), expected[i]  )