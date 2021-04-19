"""
partition.py
"""
from copy import copy
import numpy as np
from dataclasses import dataclass, field
from typing import List
from pydantic import validator, BaseModel

from .scf import scf
from .vp_nuclear import vp_nuclear
from .vp_kinetic import vp_kinetic
from .vp_hxc import vp_hxc
from .vp_overlap import vp_overlap
from .vp_surprise import vp_surprise
from .EnsCorHar import EnsCorHar
from .energy import energy
from .partition_energy import partition_energy
from .ep_nuclear import ep_nuclear
from .ep_kinetic import ep_kinetic
from .ep_hxc import ep_hxc
from .initialguessinvert import initialguessinvert
from .partition_potential import partition_potential
from ..common.coulomb import coulomb
from ..libxc.libxc import Libxc
from ..hartree.hartree import Hartree
from ..kohnsham.kohnsham import Kohnsham, KohnShamOptions

@dataclass
class V:    
    pass

@dataclass
class E:
    # Ea      : float = 0.0
    # Eb      : float = 0.0 
    # Ef      : float = 0.0
    # Tsf     : float = 0.0
    # Eksf    : List[float] = field(default_factory=list)
    # Enucf   : float = 0.0
    # Excf    : float = 0.0
    # Ecf     : float = 0.0
    # Ehf     : float = 0.0
    # Vhxcf   : float = 0.0
    # Ep_pot  : float = 0.0
    # Ep_kin  : float = 0.0
    # Ep_hxc  : float = 0.0
    # Et      : float = 0.0
    # Vnn     : float = 0.0
    # E       : float = 0.0
    # evals_a : List[float] = field(default_factory=list)
    # evals_b : List[float] = field(default_factory=list)
    # Ep_h    : float = 0.0
    # Ep_x    : float = 0.0
    # Ep_c    : float = 0.0
    pass

class PartitionOptions(KohnShamOptions):
    vp_calc_type      : str = 'component'
    hxc_part_type     : str = 'exact'
    kinetic_part_type : str = 'vonweiz'
    k_family          : str = 'gga'
    ke_func_id        : int = 5
    ke_param          : dict = {}
    ab_sym            : bool = False
    fractonal         : bool = False
    ens_spin_sym      : bool = False
    isolated          : bool = False
    fixed_q           : bool = False

    @validator('vp_calc_type')
    def vp_calc_type_values(cls, v):
        values = ['component', 'potential_inversion']
        if v not in values:
            raise ValueError(f"'vp_calc_type' must be one of the options: {values}")
        return v

    @validator('hxc_part_type')
    def hxc_part_type_values(cls, v):
        values = ['exact', 'overlap_hxc_2', 'overlap_hxc', 'overlap_xc', 'surprisal', 'hartree']
        if v not in values:
            raise ValueError(f"'hxc_part_type' must be one of the options: {values}")
        return v

    @validator('kinetic_part_type')
    def kinetic_part_type_values(cls, v):
        values = ['vonweiz','inversion','libxcke','paramke','none','twoorbital','fixed']
        if v not in values:
            raise ValueError(f"'kinetic_part_type' must be one of the options: {values}")
        return v

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
        Information about most recent inversion
    Alpha, Beta:
        Convergence Parameters
    """


    def __init__(self, grid,
                       Za, Zb,
                       pol, 
                       Nmo_a, N_a, nu_a, 
                       Nmo_b, N_b, nu_b, 
                       optPartition={}):

        #Validate options
        optPartition =  {k.lower(): v for k, v in optPartition.items()}
        for i in optPartition.keys():
            if i not in PartitionOptions().dict().keys():
                raise ValueError(f"{i} is not a valid option for KohnSham")
        optPartition = PartitionOptions(**optPartition)
        self.optPartition = optPartition

        #Grid
        self.grid = grid

        #Libxc function for fragment calculations
        # self.exchange = None
        # self.correlation = None
        # self.kinetic = None
        self.inverter = None
        self.hartree = None

        #Polarization
        self.pol = pol

        #Fragmentation 
        self.Nf = None

        #Component molecular potentials and total energies
        self.V = V()
        self.E = E()

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
        if optPartition.ab_sym and self.Za != self.Zb:
            raise ValueError("optPartition.ab_sym is set but nuclear charges are not symmetric")

        self.inversion_info = None

        #Conversion parameters
        self.Alpha = None
        self.Beta = None

        if optPartition.interaction_type == "dft":
            self.exchange = Libxc(self.grid, optPartition.xc_family, optPartition.xfunc_id)
            self.correlation = Libxc(self.grid, optPartition.xc_family, optPartition.cfunc_id)
            self.hartree = Hartree(grid,
                                    #optPartition,
                                    )
        else:
            self.exchange = 0.0
            self.correlation = 0.0
            self.hartree = 0.0

        optKS = dict( (k, getattr(optPartition, k)) for k in ('interaction_type', 
                                                            'sym', 
                                                            'fractional', 
                                                            'xfunc_id', 
                                                            'cfunc_id', 
                                                            'xc_family') if hasattr(optPartition, k) )
                                                
        #Set up kohn sham objects
        self.KSa = Kohnsham(self.grid, self.Za, 0, self.pol, self.Nmo_a, self.N_a, optKS)
        self.KSb = Kohnsham(self.grid, 0, self.Zb, self.pol, self.Nmo_b, self.N_b, optKS)

        #Figure out scale factors
        self.calc_scale_factors()

        if optPartition.kinetic_part_type == "libxcke":
            self.kinetic = Libxc(self.grid, optPartition.k_family, optPartition.ke_func_id)
        elif optPartition.kinetic_part_type == "paramke":
            self.kinetic = Paramke(self.grid, optPartition.k_family, 
                                              optPartition.ke_func_id, 
                                              optPartition.ke_param)
        
        self.calc_nuclear_potential()

#-->Methods
    def calc_scale_factors(self):
        """
        Calculates scale factors
        """
        print("Warning: If len(KS) > 1 Has not been migrated from matlab")

        self.KSa.scale = self.nu_a
        self.KSb.scale = self.nu_b

        #IF ENS_SPIN_SYM is set, then each scale factor is Reduced by a factor of 
        #two because it will be combined with an ensemble component with 
        #flipped spins
        if self.optPartition.ens_spin_sym is True:
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
            self.V.vext[:,1] = self.va + self.vb

    def mirrorAB(self):
        "Mirror fragment A to get B"

        #Mirror densities and Q functions
        self.KSb.n = self.grid.mirror(self.KSa.n).copy()
        self.KSb.Q = self.grid.mirror(self.KSa.Q).copy()

        #Energies don't need mirrored, just transfered
        self.KSb.E = copy(self.KSa.E)
        self.KSb.u = copy(self.KSa.u)

        self.KSb.veff = self.grid.mirror(self.KSa.veff).copy()
        self.KSb.vext = self.grid.mirror(self.KSa.vext).copy()

        attributes = ["vx", "vc", "vh", "vp", "ex", "ec", "eh"]
        #Mirror all potentials
        for i in attributes:
            if hasattr( self.KSa.V, i ) is True:
                setattr( self.KSb.V, i, self.grid.mirror(getattr(self.KSa.V, i)).copy() )

    def calc_protomolecule(self):
        """
        Calculate protomolecular density
        """
        
        #Evaluate sum of fragment densities and weighing functions
        self.na_frac  = np.zeros_like(self.KSa.n)
        self.nb_frac  = np.zeros_like(self.KSb.n) 
        self.na_frac += self.KSa.n * self.KSa.scale 
        self.nb_frac += self.KSb.n * self.KSb.scale

        #If we have spin symmetry in the ensemble then add
        #each density with spin flipped version
        if self.optPartition.ens_spin_sym is True:
            self.na_frac += self.grid.spinflip(self.na_frac)
            self.nb_frac += self.grid.spinflip(self.nb_frac)

        #Nf is the sum of the fragment densities
        self.nf = self.na_frac + self.nb_frac

    def calc_Q(self):
        """
        Calculate Q functions
        """ 
        np.seterr(divide='ignore', invalid='ignore')

        self.KSa.Q = self.KSa.scale * self.KSa.n / self.nf
        self.KSb.Q = self.KSb.scale * self.KSb.n / self.nf
        self.KSa.Q = np.nan_to_num(self.KSa.Q, nan=0.0, posinf=0.0, neginf=0.0)
        self.KSb.Q = np.nan_to_num(self.KSb.Q, nan=0.0, posinf=0.0, neginf=0.0)

    def vp_nuclear(self):
        vp_nuclear(self)

    def vp_kinetic(self):
        vp_kinetic(self)

    def vp_hxc(self):
        vp_hxc(self)

    def vp_overlap(self):
        vp_overlap(self)
    
    def vp_surprise(self):
        vp_surprise(self)

    def energy(self):
        energy(self)

    def partition_energy(self):
        partition_energy(self)

    def ep_nuclear(self):
        ep_nuclear(self)

    def ep_kinetic(self):
        ep_kinetic(self)

    def ep_hxc(self):
        ep_hxc(self)

    def EnsCorHar(self):
        EnsCorHar(self)

    # def get_ts_WFI(self):
    #     ts = get_ts_WFI(self)

    def partition_potential(self):
        vp = partition_potential(self)
        return vp

    def initialguessinvert(self, ispin):
        phi0, e0, vs0 = initialguessinvert(self, ispin)
        return phi0, e0, vs0

    def scf(self, optSCF={}):
        scf(self, optSCF)



