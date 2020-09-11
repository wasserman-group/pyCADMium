"""
scf.py
"""

def scf(self, optKS):
    """
    SCF method to handle self consistent field calculations
    """

    Tolerance = optKS["Tolerance"] if "Tolerance" in optKS.keys() else 10**-7
    MaxIter = optKS["MaxIter"] if "MaxIter" in optKS.keys() else 50
    Alpha = optKS["Alpha"] if "Alpha" in optKS.keys() else [0.82]
    Beta = optKS["Beta"] if "Beta" in optKS.keys() else [0.82]
    CalcType = optKS["CalcType"] if "CalcType" in optKS.keys() else "PDFT"
    Verbose = optKS["Verbose"] if "Verbose" in optKS.keys() else False
    Continue = optKS["Continue"] if "Continue" in optKS.keys() else False
    Iterative = optKS["Iterative"] if "Iterative" in optKS.keys() else True
    FreezeA = optKS["FreezeA"] if "FreezeA" in optKS.keys() else False
    FreezeB = optKS["FreezeB"] if "FreezeB" in optKS.keys() else False
    AutoTol = optKS["AutoTol"] if "Autotol" in optKS.keys() else False
    AutoTolIter = optKS["AutoTolIter"] if "AutoTolIter" in optKS.keys() else 3

    self.KSa.V["frozen"] = FreezeA   
    self.KSb.V["frozen"] = FreezeB

    if len(Alpha) == 1:
        Alpha = [Alpha, Alpha]
    elif len(Alpha) == 2:
        Alpha = Alpha
    else:
        raise ValueError ("Max length of Alpha is 2")

    self.KSa.Alpha = Alpha[0]
    self.KSb.Alpha = Alpha[1]

    if Verbose is True:
        if self.kinetic_part_type == "inversion":    
            print(f"                       Total Energy                      Inversion  \n")
            print(f"iter             A              B                  iters         optimality        res \n")
            print("------------------------------------------------------------------------------------------  \n")

        else:
            print(f"                  Total Energy            \n")
            print(f"iter              A            B              res     \n")
            print("---------------------------------------------------------\n")


    #Initial Guess Calculations
    if self.AB_SYM is True:
        #Only do calculations for fragment a
        KSab = self.KS




