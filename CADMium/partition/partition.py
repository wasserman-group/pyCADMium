"""
partition.py
"""

from pydantic import Field, constr, validator, BaseModel


class Partition():
    """ 
    Handles calculation of all efective partiions. 
    Provides an interface with libxc, hartree, and coulomb classes and functions

    Attributes
    ----------

    grid: psgrid
        PsGrid object
    interaction type: str
        {"dft", "non-interacting"}
    vp_calc_type: str
        Method for calculation of the partition potential
    xc_family: str
        Functional family
    x_func_id:  int
        Exchange functional id
    c_func_id: int 
        Correlation functional id 
    exchange:
        Lambda function for exchange functional
    correlation:
        Lambda function for correlation functional 
    k_family: str
        Kinetic energy functional family
    ke_func_id: str
        Kinetic energy functional id
    ke_param: dict
        Kinetic energy parameters
    kinetic:
        Function for non-additive kinetic energy functiona
    inverter:
        Inverter Lambda function for kinetic energy
    hartree:
        Hatree Lambda function
    hxc_part_type: str
        {"DFT", "non-interacting", "overlap approximation"}
    kinetic_part_type:
        {"von-weizacker", "inversion", "approx kinetic energy functional"}
    polarization:
        Polarization of fragments
    Nf: int
        Number of fragments
    V:
        Compotent molecular potential
    E:
        Total Energies 
    KSa, KSb: Kohn Sham Object
        Kohn Sham objects for fragments A and B
    Za, Zb: Integer
        Fragment nuclear charges
    va, vb: np.array
        Fragment potentials
    na_frac, nb_frac:
        Fragment ensembles
    nu_a, nu_b: float
        Mixing rations
    nf:
        Sum of fragment ensembles
    AB_SYM: bool
        Use of AB symmetry for homonuclear diatomics
    ENS_SPIN_SYM: bool
        Are the ensembles spin symmetric
    ISOLATED: bool
        Prescence of partition potential. Checks isolated energies
    FIXEDQ: bool
        Fixed local Q approximation
    FRACTIONAL: bool
        Allow fractional occupation of the HOMO
    inversion_info:
        Information about most recent invesrion
    Alpha, Beta:
        Convergence Parameters
    """


    grid : constr(regex="CADMium.psgrid.psgrid.Psgrid") = Field(
        "CADMium.psgrid.psgrid.Psgrid",
        description=("CADMium Prolate Spheroidal Object"),
    )

    def __init__(self, grid,
                       Za, Zb,
                       pol, 
                       Nmo_a, N_a, nu_a, Nmo_b, N_b, nu_b, 
                       optPartition,

                       interaction_type  = 'dft',
                       vp_calc_type      = 'component',
                       hxc_part_type     = 'exact',
                       kinetic_part_type = 'vonweiz', 

                       AB_SYM = False,
                       FRACTIONAL = False, 
                       ENS_SPIN_SYM = False, 
                       ISOLATED = False, 
                       FIXEDQ = False,

                       x_func_id = 1,
                       c_func_id = 12, 
                       xc_family = 'lda',
                       k_family = 'gga',
                       ke_func_id = '5',
                       ke_param = [], 
                       ):

        #Calculation Options Validators:
        if interaction_type not in ["dft", "ni"]:
            raise ValueError("Only {'dft', 'ni'} are valid options for calculation")
        else:
            self.interaction_type = interaction_type
        
        if vp_calc_type not in ["component", "potential_inversion"]:
            raise ValueError("Only {'component', 'potential_inversion'} are valid options for vp")
        else:
            self.vp_calc_type = vp_calc_type

        if hxc_part_type not in ["exact", "overlap_hxc", "overpal_xc", "surprisal", "hartree"]:
            raise ValueError("Only {'exact', 'overlap_hxc', 'overpal_xc', 'surprisal', 'hartree'}  \
                                are valid options for vp")
        else:
            self.hxc_part_type = hxc_part_type

        if kinetic_part_type not in ["vonweiz", "inversion", "libxcke", "parmke", "two_orbital", "fixed"]:
            raise ValueError("Only {'vonweiz', 'inversion', 'libxcke', 'parmke', 'two_orbital', 'fixed'} \
                                are valud options for vp")
        else:
            self.kinetic_part_type = kinetic_part_type

        #Missing type assertions

        self.grid = grid

        #xc options
        self.xc_family = xc_family
        self.x_func_id = x_func_id
        self.c_func_id = c_func_id
        
        self.exchange = None
        self.correlation = None

        #Kinetic Energy
        self.k_family = k_family
        self.ke_func_id = ke_func_id
        self.ke_param = ke_param

        self.kinetic = None
        self.inverter = None
        self.hartree = None

        self.kientic_part_type = None

        #Polarization
        self.pol = pol

        #Fragmentation 
        self.Nf = None

        #Component molecular potentials and total energies
        self.V = None
        self.E = None

        #Kohn Sham objects
        self.KSa = None
        self.KSb = None

        #Fragment nuclear charges and potentials
        self.Za, self.Zb = Za, Zb
        self.va, self.vb = None, None

        #Fragment ensembles, mixing rations, sum of fragment ensembles
        self.na_frac, self.nb_fac = None, None
        self.nu_a, self.nu_b = nu_a, nu_b
        self.Nmo_a, self.Nmo_b = Nmo_a, Nmo_b
        self.nf = None

        #Flags
        AB_SYM = AM_SYM
        ENS_SPIN_SYM = ENS_SPIN_SYM
        ISOLATED = ISOLATED
        FIXEDQ = FIXEDQ
        FRACTIONAL = FRACTIONAL

        inversion_info = False

        #Conversion parameters
        Alpha = None
        Beta = None
    