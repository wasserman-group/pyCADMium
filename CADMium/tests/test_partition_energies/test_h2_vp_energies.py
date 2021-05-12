import numpy as np
from CADMium import Pssolver, Psgrid, Partition, Inverter


def test_partition_energies():

    a = 1.466/2
    #Nuclear charge for fragments A and B
    Za, Zb = 1,1
    #Set polarization 1-Unpolarized, 2-Polarized|
    pol = 2
    #Fragment a electrons [alpha, beta]
    Nmo_a = [[1,0]] #Number of molecular orbitals to calculate
    N_a   = [[1,0]]
    #Ensemble mix
    nu_a = 1
    #Fragment b electrons
    Nmo_b = [[1,0]]
    N_b   = [[1,0]]
    #Ensemble mix
    nu_b = 1

    #Molecular elctron configuration
    Nmo_m = [[1,1]]
    N_m   = [[1,1]]

    #Set up grid
    NP = 4
    NM = [4,4]
    L = np.arccosh(12/a)
    loc = np.array(range(-4,5)) #Stencil outline

    grid = Psgrid(NP, NM, a, L, loc)
    grid.initialize()

    part = Partition(grid, Za, Zb, pol, Nmo_a, N_a, nu_a, Nmo_b, N_b, nu_b, {  "AB_SYM"       : True,
                                                                            "ENS_SPIN_SYM" : True,  
                                                                            "kinetic_part_type" : "inversion",
                                                                            "k_family" : "gga",
                                                                            "ke_func_id" : 500,
                                                                                })

    #Setup inverter object
    mol_solver = Pssolver(grid, Nmo_m, N_m)
    part.inverter = Inverter(grid, mol_solver, {  "AB_SYM"         : True,
                                                "ENS_SPIN_SYM"   : True,  
                                                "use_iterative"  : False,
                                                "invert_type"    : "orbitalinvert",
                                                "DISP"           : False,  
                                                })

    part.optPartition.isolated   = False

    part.scf({"disp"       : False,
            "alpha"      : [0.6],
            "max_iter"   : 200,
            "e_tol"      : 1e-9,
            "iterative"  : False,
            "continuing" : False})
    
    expected = {'Ea': -0.45648534732914625,
                'Eb': -0.45648534732914625,
                'Ef': -0.9129706946582925,
                'Tsf': 1.229793692709058,
                'Eksf': np.array([[-0.74204158,  0.        ]]),
                'Enucf': -2.195184927508663,
                'Exf': -0.5953289679034595,
                'Ecf': -0.0464215259607141,
                'Ehf': 0.6941710340054861,
                'Vhxcf': 0.5414671244122083,
                'Ep_pot': -1.3335441021724597,
                'Ep_kin': -0.15202727101015334,
                'Ep_h': 0.5801249630265373,
                'Ep_x': 0.04570980033028593,
                'Ep_c': -0.046872964090811736,
                'Ep_hxc': 0.5789617992660114,
                'Ep': -0.9066095739166016,
                'Et': -1.819580268574894,
                'Vnn': 0.6821282401091405,
                'E': -1.1374520284657534,
                'evals_a': np.array([], dtype=np.float64),
                'evals_b': np.array([], dtype=np.float64)}

    for i in part.E.__dict__:
        if i.startswith("__") is False:
            if type(getattr(part.E, i)) == np.ndarray:
                assert np.isclose(  getattr(part.E, i).all(), expected[i].all()  )
            else:
                assert np.isclose(  getattr(part.E, i), expected[i]  )