"""
partition.py
"""
import numpy as np


from .scf import scf
from ..common.coulomb import coulomb
from ..libxc.libxc import Libxc
from ..hartree.hartree import Hartree
from ..kohnsham.kohnsham import Kohnsham

class V:
    def __iter__(self):
        for attr, value in self.__dict__.iteritems():
            yield attr, value
            
    pass

class E:
    pass

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


    # grid : constr(regex="CADMium.psgrid.psgrid.Psgrid") = Field(
    #     "CADMium.psgrid.psgrid.Psgrid",
    #     description=("CADMium Prolate Spheroidal Object"),
    # )

    def __init__(self, grid,
                       Za, Zb,
                       pol, 
                       Nmo_a, N_a, nu_a, 
                       Nmo_b, N_b, nu_b, 
                       optPartition={}):

        optPartition["interaction_type"] = optPartition["interaction_type"] if "interaction_type" in optPartition.keys() else "dft"
        optPartition["vp_calc_type"] = optPartition["vp_calc_type"] if "vp_calc_type" in optPartition.keys() else "component"
        optPartition["hxc_part_type"] = optPartition["hxc_part_type"] if "hxc_part_type" in optPartition.keys() else "exact"
        optPartition["kinetic_part_type"] = optPartition["kinetic_part_type"] if "kinetic_part_type" in optPartition.keys() else "vonweiz"

        optPartition["AB_SYM"] = optPartition["AB_SYM"] if "AB_SYM" in optPartition.keys() else False
        optPartition["FRACTIONAL"] = optPartition["FRACTIONAL"] if "FRACTIONAL" in optPartition.keys() else False
        optPartition["ENS_SPIN_SYM"] = optPartition["ENS_SPIN_SYM"] if "ENS_SPIN_SYM" in optPartition.keys() else False
        optPartition["ISOLATED"] = optPartition["ISOLATED"] if "ISOLATED" in optPartition.keys() else False
        optPartition["FIXEDQ"] = optPartition["FIXEDQ"] if "FIXEDQ" in optPartition.keys() else False
        
        optPartition["x_func_id"] = optPartition["x_func_id"] if "x_func_id" in optPartition.keys() else 1
        optPartition["c_func_id"] = optPartition["c_func_id"] if "c_func_id" in optPartition.keys() else 12
        optPartition["xc_family"] = optPartition["xc_family"] if "xc_family" in optPartition.keys() else "lda"
        optPartition["k_family"] = optPartition["k_family"] if "k_family" in optPartition.keys() else "gga"
        optPartition["ke_func_id"] = optPartition["ke_func_id"] if "ke_func_id" in optPartition.keys() else 5
        optPartition["ke_param"] = optPartition["ke_param"] if "ke_param" in optPartition.keys() else {}
        
        self.optPartition = optPartition

        #Calculation Options Validators:
        if self.optPartition["interaction_type"]  not in ["dft", "ni"]:
            raise ValueError("Only {'dft', 'ni'} are valid options for calculation")
        
        if self.optPartition["vp_calc_type"]  not in ["component", "potential_inversion"]:
            raise ValueError("Only {'component', 'potential_inversion'} are valid options for vp")

        if self.optPartition["hxc_part_type"] not in ["exact", "overlap_hxc", "overpal_xc", "surprisal", "hartree"]:
            raise ValueError("Only {'exact', 'overlap_hxc', 'overpal_xc', 'surprisal', 'hartree'}  \
                                are valid options for vp")

        if self.optPartition["kinetic_part_type"] not in ["vonweiz", "inversion", "libxcke", "parmke", "two_orbital", "fixed"]:
            raise ValueError("Only {'vonweiz', 'inversion', 'libxcke', 'parmke', 'two_orbital', 'fixed'} \
                                are valud options for vp")

        #Missing type assertions

        self.grid = grid

        #xc options
        # self.xc_family = xc_family
        # self.x_func_id = x_func_id
        # self.c_func_id = c_func_id

        #Libxc function for fragment calculations
        self.exchange = None
        self.correlation = None
        
        #Kinetic Energy
        # self.k_family = k_family
        # self.ke_func_id = ke_func_id
        # self.ke_param = ke_param

        self.kinetic = None
        self.inverter = None
        self.hartree = None

        self.kientic_part_type = None

        #Polarization
        self.pol = pol

        #Fragmentation 
        self.Nf = None

        #Component molecular potentials and total energies
        self.V = V
        self.E = E

        #Kohn Sham objects
        self.KSa = None
        self.KSb = None

        #Fragment nuclear charges and potentials
        self.Za, self.Zb = Za, Zb
        self.va, self.vb = None, None

        #Fragment ensembles, mixing rations, sum of fragment ensembles
        self.na_frac, self.nb_fac = None, None
        self.nu_a, self.nu_b = nu_a, nu_b
        self.N_a = np.array(N_a)
        self.N_b = np.array(N_b)
        self.Nmo_a = np.array(Nmo_a)
        self.Nmo_b = np.array(Nmo_b)
        self.nf = None

        #Sanity Check
        if self.optPartition["AB_SYM"] and self.Za != self.Zb:
            raise ValueError("AB_SYM is set but nuclear charges are not symmetric")

        inversion_info = False

        #Conversion parameters
        self.Alpha = None
        self.Beta = None

        if self.optPartition["interaction_type"] == "dft":
            self.exchange = Libxc(self.grid, self.optPartition["xc_family"], self.optPartition["x_func_id"])
            self.correlation = Libxc(self.grid, self.optPartition["xc_family"], self.optPartition["c_func_id"])
            self.hartree = Hartree(grid,
                                    #optPartition,
                                    )
        else:
            self.exchange = 0.0
            self.correlation = 0.0
            self.hartree = 0.0

    
        #Set up kohn sham objects
        self.KSa = Kohnsham(self.grid, self.Za, 0, self.pol, self.Nmo_a, self.N_a, optPartition)
        self.KSb = Kohnsham(self.grid, 0, self.Zb, self.pol, self.Nmo_b, self.N_b, optPartition)

        #Figure out scale factors
        self.calc_scale_factors()

        if self.optPartition["kinetic_part_type"] == "libxcke":
            self.kinetic = Libxc(self.grid, self.k_family, self.ke_func_id)
        elif self.optPartition["kinetic_part_type"] == "paramke":
            self.kinetic = Paramke(self.grid, self.k_family, self.ke_func_id, self.ke_param)
        
        self.calc_nuclear_potential()

    #Methods
    def calc_scale_factors(self):
        """
        Calculates scale factors
        """
        #print("Warning: If len(KS) > 1 Has not been migrated from matlab")

        self.KSa.scale = self.nu_a
        self.KSb.scale = self.nu_b

        #IF ENS_SPIN_SYM is set, then each scale factor is Reduced by a factor of 
        #two because it will be combined with an ensemble component with 
        #flipped spins

        if self.optPartition["ENS_SPIN_SYM"] is True:
            self.KSa.scale = self.KSa.scale / 2.0
            self.KSb.scale = self.KSb.scale / 2.0

    def calc_nuclear_potential(self):
        """
        Calculate external nuclear potentials
        """

        self.va = coulomb(self.grid, self.Za, 0)
        self.vb = coulomb(self.grid, 0, self.Zb)

        self.V.vext = np.zeros((self.va.shape[0], self.pol))

        self.V.vext[:,0] = self.va + self.vb
        if self.pol == 2:   
            self.V.vext[0,1] = self.va + self.vb

    def mirrorAB(self):
        "Mirror fragment A to get B"

        #Mirror densities and Q functions
        self.KSb.n = self.grid.mirror(self.KSa.n)
        self.KSb.Q = self.grid.mirror(self.KSa.Q)

        #Energies don't need mirrored, just transfered
        self.KSb.E = self.KSa.E
        self.KSb.u = self.KSa.u

        self.KSb.veff = self.grid.mirror(self.KSa.veff)
        self.KSb.vext = self.grid.mirror(self.KSa.vext)

        #Mirror all the potentials
        for attribute in self.KSa.V.__dict__.keys():
            if "__" not in attribute:
                self.KSb.V.attribute = self.grid.mirror(self.KSa.V.attribute)

    def scf(self):
        scf(self, optPartition=self.optPartition)



