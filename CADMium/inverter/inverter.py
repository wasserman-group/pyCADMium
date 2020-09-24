"""
inverter.py
"""

class Inverter():

    def __init__(self, grid, solver, optInversion):
        """
        Inverter handles inversion algorithms

        Parameters 
        ----------

        Attributes
        ----------
        vs : numpy array
            Kohn Sham Potential
        us : float
            Chemical Potential
        ts_WFI :
            Kinetic energy density using laplacian
        ts_WFII :
            Kinetic energy density using gradient
        Ts : float
            Kinetic Energy


        Defaults
        --------
        Convergence and Algorithm Parameters

        Use_Iterative : logical
            Set to True to use iterative sovler with orbital invert or eigesolveinvert.
            Otherwise direct solver will be used which may take a long time for large
            grids.

        Kval : floag
            Parameter for iterative preconditioner

        TolLinsolve : float
            Tolerance for iterative solution to linear solve

        MaxIterLinsolve: int
            Maximum number of iterations for solution to linear solve

        TolInvert: float
            Requested tolerance for inversion

        MaxIterInvert: int
            Maximum iterations for inversions

        ResFactor: 1e0
            Determines choice of update in orbitalinvert
            
                    """


        optInversion["AB_SYM"] = optInversion["AB_SYM"] if "AB_SYM" in optInversion.keys() else False
        optInversion["ENS_SPIN_SYM"] = optInversion["ENS_SPIN_SYM"] if "ENS_SPIN_SYM" in optInversion.keys() else False
        optInversion["USE_ITERATIVER"] = optInversion["USE_ITERATIVE"] if "USE_ITERATIVE" in optInversion.keys() else False
        optInversion["VERBOSE"] = optInversion["VERBOSE"] if "VERBOSE" in optInversion.keys() else False
        
        optInversion["AVOIDLOOP"] = optInversion["AVOIDLOOP"] if "AVOIDLOOP" in optInversion.keys() else False

        optInversion["invert_type"] = optInversion["invert_type"] if "invert_type" in optInversion.keys() else "wuyang"
        optInversion["Kval"] = optInversion["Kval"] if "Kval" in optInversion.keys() else -1e-12
        optInversion["TolLinsolve"] = optInversion["TolLinsolve"] if "TolLinsolve" in optInversion.keys() else 1e-2
        optInversion["TolInvert"] = optInversion["TolInvert"] if "TolInvert" in optInversion.keys() else 1e-12
        optInversion["MaxIterLinsolve"] = optInversion["MaxIterLinsolve"] if "MaxIterLinsolve" in optInversion.keys() else 2000
        optInversion["MaxIterInvert"] = optInversion["MaxIterInvert"] if "MaxIterInvert" in optInversion.keys() else 20
        optInversion["ResFactor"] = optInversion["ResFactor"] if "ResFactor" in optInversion.keys() else 1e0

        self.optInversion = optInversion

        self.grid = grid
        self.solver = solver

        self.vs = None
        self.us = None
        self.ts_WFI = None
        self.ts_WFII = None

    def invert(self, n0, vs0, phi0=[], e0=[], ispin=[], Qi=[]):
        """
        Does the inverstion
        """

        if self.optInversion["invert_type"] == "wuwang":
            flag, output = self.linresponse(n0, vs0, ispin)

        elif self.optInversion["invert_type"] == "simple":
            self.simple(n0, vs0, ispin)

        elif self.optInversion["invert_type"] == "orbitalinvert":
            flag, output = self.orbitalinvert(n0, vs0, phi0, e0, ispin)

        elif self.optInversion["invert_type"] == "qinvert":
            flag, output = self.qinvert(n0, vs0, phi0, e0, ispin, Qi)

        elif self.optInversion["invert_type"] == "eigensolveinvert":
            flag, output = self.eigensolveinvert(n0, vs0, ispin)

        elif self.optInversion["invert_type"] == "test":
            flag, output = self.test(n0, vs0, phi0, e0, ispin)


