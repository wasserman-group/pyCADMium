"""
scf.py
"""

import numpy as np

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


    self.KSa.frozen = optPartition["FreezeA"]   
    self.KSb.frozen = optPartition["FreezeB"]

    if len(optPartition["Alpha"]) == 1:
        self.Alpha = [optPartition["Alpha"], optPartition["Alpha"]]
    elif len(optPartition["Alpha"]) == 2:
        self.Alpha = optPartition["Alpha"]
    else:
        raise ValueError ("Max length of Alpha is 2")

    self.KSa.Alpha = self.Alpha[0]
    self.KSb.Alpha = self.Alpha[1]

    if optPartition["Verbose"] is True:
        if self.kinetic_part_type == "inversion":    
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




