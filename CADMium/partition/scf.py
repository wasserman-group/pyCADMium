"""
scf.py
"""

import numpy as np
import sys

def scf(self, optPartition):
    """
    SCF method to handle self consistent field calculations
    """

    optPartition["Tolerance"] = optPartition["Tolerance"] if "Tolerance" in optPartition.keys() else 10**-7
    optPartition["MaxIter"] = optPartition["MaxIter"] if "MaxIter" in optPartition.keys() else 50
    optPartition["Alpha"] = optPartition["Alpha"] if "Alpha" in optPartition.keys() else [0.82]
    optPartition["Beta"] = optPartition["Beta"] if "Beta" in optPartition.keys() else [0.82]

    optPartition["CalcType"] = optPartition["CalcType"] if "CalcType" in optPartition.keys() else "PDFT"
    optPartition["Verbose"] = optPartition["Verbose"] if "Verbose" in optPartition.keys() else False
    optPartition["CONTINUE"] = optPartition["CONTINUE"] if "CONTINUE" in optPartition.keys() else False
    optPartition["ITERATIVE"] = optPartition["ITERATIVE"] if "ITERATIVE" in optPartition.keys() else True
    optPartition["AutoTol"] = optPartition["AutoTol"] if "Autotol" in optPartition.keys() else False
    optPartition["AutoTolIter"] = optPartition["AutoTolIter"] if "AutoTolIter" in optPartition.keys() else 3

    optPartition["FreezeA"] = optPartition["FreezeA"] if "FreezeA" in optPartition.keys() else False
    optPartition["FreezeB"]  = optPartition["FreezeB"] if "FreezeB" in optPartition.keys() else False

    self.optPartition = optPartition


    self.KSa.V.frozen = optPartition["FreezeA"]   
    self.KSb.V.frozen = optPartition["FreezeB"]

    if len(optPartition["Alpha"]) == 1:
        self.Alpha = [optPartition["Alpha"][0], optPartition["Alpha"][0]]
    elif len(optPartition["Alpha"]) == 2:
        self.Alpha = optPartition["Alpha"]
    else:
        raise ValueError ("Max length of Alpha is 2")

    self.KSa.Alpha = self.Alpha[0]
    self.KSb.Alpha = self.Alpha[1]

    if optPartition["Verbose"] is True:
        if self.optPartition["kinetic_part_type"] == "inversion":    
            print(f"                       Total Energy                      Inversion  \n")
            print(f"iter             A              B                  iters         optimality        res \n")
            print("------------------------------------------------------------------------------------------  \n")

        else:
            print(f"                  Total Energy            \n")
            print(f"iter              A            B              res     \n")
            print("---------------------------------------------------------\n")


    #Initial Guess Calculations
    if optPartition["AB_SYM"] is True:
        #Only do calculations for fragment a
        KSab = [self.KSa]

    else:
        KSab = [self.KSa, self.KSb]


    for i_KS in KSab:
        if optPartition["CONTINUE"] is True:
            #If we are continuiing a calculation, check that we have an input density
            assert i_KS.n != None, "CONTINUE flag set but there is no input density"
        
        else:
            #We need initial guesses
            i_KS.vext = np.zeros_like(i_KS.vnuc)
            i_KS.vhxc = np.zeros_like(i_KS.vnuc)
            i_KS.set_effective_potential()

            #Initial guess for effective potential is just nuclear potential
            nout = i_KS.calc_density(optPartition["ITERATIVE"])
            i_KS.n = nout
            i_KS.calc_chempot()

    if optPartition["AB_SYM"] == True:
        self.mirrorAB()

    #Form protomolecule and calculate Q-functions
    self.calc_protomolecule()
    self.calc_Q()

    #Start up the scf loop
    dif        = 10.0
    old_E      = 0.0
    old_nf     = self.nf
    vp         = np.zeros_like((self.grid.Nelem, self.pol))
    iterations = 1
    failures   = 0
    STOP       = False

    if optPartition["AutoTol"] is True:
        min_dif = 10.0
        num_iter_not_min = 0
    
    while (dif > optPartition["Tolerance"] 
           or (optPartition["AutoTol"] is True and num_iter_not_min < optPartition["AutoTolIter"])) \
           and iterations <= optPartition["MaxIter"]                                                \
           and STOP is False:

        #Set flags to avoid dead loops in inversion
        if optPartition["kinetic_part_type"] == "inversion" and optPartition["ISOLATED"] is False:
            if iterations == 1:
                self.inverter.optInversion["AVOIDLOOP"] = True
            else:
                self.inverter.optInversion["AVOIDLOOP"] = False
            
        for iKS in KSab:
            #Calculate new local vhxc potentials
            vhxc_old = iKS.vhxc
            iKS.calc_hxc_potential()
            iKS.vhxc = 1.0 * iKS.vhxc + 0.0 * vhxc_old
        
        if self.optPartition["AB_SYM"]:
            self.mirrorAB()

        if not self.optPartition["ISOLATED"]:
            #Calculate the partition potential
            vp = 0.0 * vp + 1.0 * self.partition_potential()

        else: 
            #Otherwise fill in zeros
            vp = np.zeros_like(self.nf)
            for i_KS in KSab:
                i_KS.V.vp = np.zeros_like(self.nf)

        for i_KS in KSab:
            if i_KS.V.frozen is not True:
                #Each fragment feels the effective potential
                #as an external potential
                if optPartition["CalcType"] == "pdft":
                    #Global partition optential
                    i_KS.vext = vp
                elif optPartition["CalcType"] == "sdft":
                    #Sybsystem-DFT
                    #Fragment embedding poetntials
                    i_KS.vext = i_KS.V.vp

            #Set new effective potential
            i_KS.set_effective_potential()

            #Calculate new density
            nout = i_KS.calc_density(self.optPartition["ITERATIVE"], dif)

            #Get the new chemical potential
            i_KS.calc_chempot()

            #Set new density with linear mixing
            i_KS.n = (1-i_KS.Alpha) * i_KS.n + i_KS.Alpha * nout

            #Calculate each fragment energy
            i_KS.energy()

        if optPartition["AB_SYM"]  is True:
            self.mirrorAB()

        #Form promolecule and calculate q functions
        self.calc_protomolecule()
        self.calc_Q()
        self.energy()

        #Convergence check
        dif_E  = np.abs( np.sum(self.E.E - old_E) / self.E.E )
        dif_nf = self.grid.integrate(np.abs(self.nf-old_nf)) / self.grid.integrate(np.abs(self.nf))
        dif    = max(dif_E, dif_nf)
        old_E  = self.E.E
        old_nf = self.nf

        if self.optPartition["AutoTol"] is True:
            if dif < min_dif:
                num_iter_not_min = 0
                min_dif = dif

            else:
                num_iter_not_min = num_iter_not_min + 1

        if optPartition["kinetic_part_type"] == "inversion" or \
            optPartition["vp_calc_type"] == "vp_calc_type" and self.optPartition["ISOLATED"]:
        
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

        #Display
        if self.optPartition["Verbose"] is True:
            if self.optPartition["kinetic_part_type"] == "inversion" or \
                self.optPartition["vp_calc_type"] == "potential inversion" and self.optPartition["ISOLATED"]:

                print(f" {iterations}  {self.E.Ea}  {self.E.Eb}  {max_nfev}  {max_optimality} {dif}")


            else:
                print(f" {iterations}  {self.E.Ea}  {self.E.Eb}   {dif} ")



    if dif < self.optPartition["Tolerance"]:
        flag = True

    else:
        flag = False

    for i_KS in KSab:
        delattr(i_KS.V, "frozen") 


    for i in self.E.__dict__:
        print(i, getattr(self.E, str(i)))


        





