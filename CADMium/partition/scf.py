"""
scf.py
"""

import numpy as np
from warnings import warn
from pydantic import validator, BaseModel
from typing import List

class PartitionSCFOptions(BaseModel):
    e_tol               : float = 1e-7
    max_iter            : int = 50
    alpha               : List[float] = [0.82]
    beta                : List[float] = [0.82]
    calc_type           : str = 'pdft'
    disp                : bool = False
    continuing          : bool = False
    from_target_density : bool = False
    iterative           : bool = True
    isolated            : bool = False
    auto_tol            : bool = False
    auto_tol_iter       : int = 3
    freeze_a            : bool = False
    freeze_b            : bool = False 

    @validator('calc_type') 
    def calc_type_values(cls, v):
        values = ['pdft', 'sdft']
        if v not in values:
            raise ValueError(f"'calc_type' must be one of the options: {values}")
        return v


def scf(self, optSCF={}):
    """
    SCF method to handle self consistent field calculations
    """

    #---> Validate options
    optSCF =  {k.lower(): v for k, v in optSCF.items()}
    for i in optSCF.keys():
        if i not in PartitionSCFOptions().dict().keys():
            raise ValueError(f"{i} is not a valid option for KohnSham")
    optSCF = PartitionSCFOptions(**optSCF)
    optSCF.isolated = self.optPartition.isolated

    #---> Sanity Check on options
    if self.inverter != None:
        if self.optPartition.ab_sym != self.inverter.optInv.ab_sym:
            warn("Warning optPartition.ab_sym != opt.Inverter.ab_sym. Is this the inteded behaviour")
        if self.optPartition.ens_spin_sym != self.inverter.optInv.ens_spin_sym:
            warn("Warning optPartition.ens_spin_sym != opt.Inverter.ab_sym. Is this the inteded behaviour")

    #---> Do we freeze fragments?
    self.KSa.V.frozen = optSCF.freeze_a   
    self.KSb.V.frozen = optSCF.freeze_b
    if self.ens:
        self.KSA.V.frozen = optSCF.freeze_a   
        self.KSB.V.frozen = optSCF.freeze_b

    #---> SCF Mix
    if len(optSCF.alpha) == 1:
        self.Alpha = [optSCF.alpha[0], optSCF.alpha[0]]
    elif len(optSCF.alpha) == 2:
        self.Alpha = optSCF.alpha
    else:
        raise ValueError ("Max length of Alpha is 2")
    self.KSa.Alpha = self.Alpha[0]
    self.KSb.Alpha = self.Alpha[1]
    if self.ens:
        self.KSA.Alpha = self.Alpha[0]
        self.KSB.Alpha = self.Alpha[1]

    if optSCF.disp is True:
        if self.optPartition.isolated:
            print("----> Begin SCF calculation for *Isolated* Fragments\n")
        else:
            print("----> Begin SCF calculation for *Interacting* Fragments\n")

        if self.optPartition.kinetic_part_type == "inversion":    
            print(f"                Total Energy (a.u.)                                Inversion                \n")
            print(f"                __________________                ____________________________________     \n")
            print(f"Iteration         A              B                  iters      optimality        res       \n")
            print("___________________________________________________________________________________________ \n")

        else:
            print(f"                Total Energy (a.u.)       \n")
            print(f"                __________________        \n")
            print(f"Iteration         A            B              res     \n")
            print("_______________________________________________________\n")


    #---> Initial Guess
    if self.optPartition.ab_sym is True:
        #Only do calculations for fragment a
        if not self.ens:
            KSab = [self.KSa]
        else:
            KSab = [self.KSa, self.KSA]

    else:
        if not self.ens:
            KSab = [self.KSa, self.KSb]
        else:
            KSab = [self.KSa, self.KSA, self.KSb, self.KSB]

    for i_KS in KSab:
        if optSCF.continuing is True:
            assert i_KS.n is not None, "optSCF.continuing is True set but there is no input density"
        
        else:
            #We need initial guesses
            i_KS.vext = np.zeros_like(i_KS.vnuc)
            i_KS.vhxc = np.zeros_like(i_KS.vnuc)
            i_KS.set_effective_potential() # Sets nuclear potential
            nout = i_KS.calc_density(optSCF.iterative)
            i_KS.n = nout
            i_KS.calc_chempot()

    if self.optPartition.ab_sym == True:
        self.mirrorAB()

    #---> Form protomolecule and calculate Q-Functions
    if not optSCF.from_target_density: 
        self.calc_protomolecule()
    self.calc_Q()

    #---> SCF Initialization
    dif        = 10.0
    old_E      = 0.0
    old_nf     = self.nf
    iterations = 1
    inversionfailures   = 0
    STOP       = False

    if optSCF.auto_tol is True:
        min_dif = 10.0
        num_iter_not_min = 0

    #-----> SCF Procedure Begins
    while (dif > optSCF.e_tol 
           or (optSCF.auto_tol is True and num_iter_not_min < self.optPartition["AutoTolIter"])) \
           and iterations <= optSCF.max_iter                                                     \
           and STOP is False:

    #-----> Sanity Checks
        if (self.optPartition.kinetic_part_type == "inversion" and optSCF.isolated == False):
            if iterations == 1:
                self.inverter.optInv.avoid_loop = True
            else:
                self.inverter.optInv.avoid_loop = False

    #-----> Calculate V_Hxc
        for iKS in KSab:
            vhxc_old = iKS.vhxc
            iKS.calc_hxc_potential()
            iKS.vhxc = 1.0 * iKS.vhxc + 0.0 * vhxc_old
        
        if self.optPartition.ab_sym:
            self.mirrorAB()

    #-----> Calculate Partition Potential
        if not optSCF.isolated:
            vp = self.partition_potential()
        else: 
            vp = np.zeros_like(self.nf)
            for i_KS in KSab:
                i_KS.V.vp = np.zeros_like(self.nf)

    #-----> Add effective potentials to fragments
            #Calculate new density
        for i_KS in KSab:
            if i_KS.V.frozen is not True:
                if optSCF.calc_type == "pdft": #Global partition optential
                    i_KS.vext = vp
                elif optSCF.calc_type == "sdft": #Sybsystem-DFT. Fragment embedding potentials
                    i_KS.vext = i_KS.V.vp

            #Set new effective potential
            i_KS.set_effective_potential()
            #Calculate new density
            nout = i_KS.calc_density(optSCF.iterative, dif)
            #Get the new chemical potential
            i_KS.calc_chempot()
            #Set new density with linear mixing
            i_KS.n = (1-i_KS.Alpha) * i_KS.n + i_KS.Alpha * nout
            #Calculate each fragment energy
            i_KS.energy()

        if self.optPartition.ab_sym  is True:
            self.mirrorAB()
    
        #Form promolecule and calculate q functions
        if not optSCF.from_target_density: 
            self.calc_protomolecule()        
        self.calc_Q()
        self.energy()

    #-----> Convergence Check
        dif_E  = np.abs( np.sum(self.E.E - old_E) / self.E.E )
        dif_nf = np.max(self.grid.integrate(np.abs(self.nf-old_nf)) / self.grid.integrate(np.abs(self.nf)))
        dif    = max(dif_E, dif_nf)
        old_E  = self.E.E
        old_nf = self.nf

        if optSCF.auto_tol is True:
            if dif < min_dif:
                num_iter_not_min = 0
                min_dif = dif
            else:
                num_iter_not_min = num_iter_not_min + 1

        if (self.optPartition.kinetic_part_type == "inversion" or \
            self.optPartition.vp_calc_type == "potential_inversion") and not optSCF.isolated: 
    
            max_nfev = 0
            max_optimality = 1e-16

            for ii in range(self.inversion_info.shape[0]):
                for jj in range(self.inversion_info.shape[1]):

                    max_nfev = max_nfev if self.inversion_info[ii,jj].nfev < max_nfev else self.inversion_info[ii,jj].nfev
                    max_optimality = max_optimality if self.inversion_info[ii,jj].nfev < max_optimality else self.inversion_info[ii,jj].optimality

                    if self.inversion_info[ii,jj].optimality > 1:
                        inversionfailures += 1
                    else:
                        inversionfailures = 0

                    if inversionfailures > 1:
                        raise SystemExit("Too many inversion failures. Stopping")

    #---> Display
        if optSCF.disp is True:
        # if True:
            if (self.optPartition.kinetic_part_type == "inversion" or \
            self.optPartition.vp_calc_type == "potential_inversion") and not optSCF.isolated: 
                print(f"  {iterations:3d}          {self.E.Ea:+10.5f}      {self.E.Eb:+10.5f}           {max_nfev:3d}       {max_optimality:+7.3e}      {dif:+7.3e}")

            else:
                print(f"  {iterations:3d}         {self.E.Ea:10.5f}   {self.E.Eb:10.5f}       {dif:7.3e} ")

        iterations += 1
        # if iterations == 2 and self.optPartition.isolated is False:
        #     sys.exit()
        if iterations == optSCF.max_iter:
            warn("SCF Warning: Max number of Iterations Surpassed. Desired convergence may have not been achieved") 

    # if dif < optSCF.e_tol:
    #     flag = True

    # else:
    #     flag = False