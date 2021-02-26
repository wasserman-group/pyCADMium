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
    NM = [4,4]
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
                                                "disp"           : False, 
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
    part.scf({"disp"       : False,
            "alpha"      : [0.3],
            "max_iter"   : 200,
            "e_tol"      : 1e-7,
            "continuing" : True, 
            "iterative"  : False})

    expected = {'Ea': -14.64239075876027,
                'Eb': -14.64239075876027,
                'Ef': -29.28478151752054,
                'Tsf': 25.34937178021555,
                'Eksf': np.array([[-16.72707747]]),
                'Enucf': -63.73131825287787,
                'Exf': -4.527858632036746,
                'Ecf': -0.4495664797297575,
                'Ehf': 14.074590066908286,
                'Vhxcf': 21.60485817523808,
                'Ep': -3.614959527515154,
                'Ep_pot': -7.1522854050404625,
                'Ep_kin': 0.0642644656425091,
                'Ep_hxc': 3.4730614118827994,
                'Et': -32.8997410450357,
                'Vnn': 3.5382574082264484,
                'E': -29.361483636809247,
                'evals_a': np.array([], dtype=np.float64),
                'evals_b': np.array([], dtype=np.float64),
                'Ep_h': 3.519875684463109,
                'Ep_x': -0.04185713118930234,
                'Ep_c': -0.004957141391007558}

    for i in part.E.__dict__:
        if i.startswith("__") is False:
            if type(getattr(part.E, i)) == np.ndarray:
                assert np.isclose(  getattr(part.E, i).all(), expected[i].all()  )
            else:
                assert np.isclose(  getattr(part.E, i), expected[i]  )