"""
vp_kinetic.py
"""

import sys
import numpy as np

def vp_kinetic(self):
    """
    Calculate Kinetic energy components of vp
    """

    #Select method for calculating kinetic energy
    if self.optPartition["kinetic_part_type"] == "vonweiz":

        #Use the von-weizacker inversion
        #May break if only one density 
        self.V.vt = (-0.5 * self.grid.elap @ (self.nf**0.5)) / (self.nf**0.5) / (self.grid.w)

        #Evaluate kinetic energy for integer ocupations
        #Densities
        #May break if only one density
        self.KSa.V.vt = (-0.5 * self.grid.elap @ (self.KSa.n**0.5)) / (self.KSa.n**0.5) / (self.grid.w)
        self.KSb.V.vt = (-0.5 * self.grid.elap @ (self.KSb.n**0.5)) / (self.KSb.n**0.5) / (self.grid.w)

        for i in range(self.KSa.V.vt):
            if np.isnan(self.KSa.V.vt[i]) is True:
                self.KSa.V.vt[i] = 0.0
            if np.isnan(self.KSb.V.vt[i]) is True:
                self.KSb.V.vt[i] = 0.0

    elif self.optPartition["kinetic_part_type"] == "libxcke" or \
         self.optPartition["kinetic_part_type"] == "paramke":
        
        #Use a kinetic energy functional from libxc
        #Evaluate molecular kinetic energy
        raise ValueError("LibxcKe nor Paramke options are not avaliable")


    elif self.optPartition["kinetic_part_type"] == "inversion":

        #Find kinetic energy fucntional derivative for fragmetns using Euler equation
        #Fragments using euler equations
        # u = max(self.KSa.u, self.KSb.u)
        homos_a = []
        homos_b = []
        u = -1 / np.spacing(1) * np.ones((1, self.pol))
        homos_a.append(u)
        homos_b.append(u)

        #Check the values of each solver's homos and pick the largest
        for i in range(self.Nmo_a.shape[0]):
            for j in range(self.Nmo_a.shape[1]):
                if self.KSa.solver[i,j].homo is not None:
                    homos_a.append(self.KSa.solver[i,j].homo)
                if self.KSb.solver[i,j].homo is not None:
                    homos_b.append(self.KSb.solver[i,j].homo)
        u = max(max(homos_a), max(homos_b))

        nspin = self.pol

        if self.optPartition["ENS_SPIN_SYM"]:
            u = max(u)
            nspin = 1

        self.KSa.V.vt = np.ones((self.grid.Nelem, 1)) * u - self.KSa.veff
        self.KSb.V.vt = np.ones((self.grid.Nelem, 1)) * u - self.KSb.veff

        for ispin in range(nspin):
            ntarget = self.nf[:, ispin]
            if self.pol == 2:
                ntarget = 2 * ntarget

        
        phi0, e0, vs0 = self.initialguessinvert(ispin)

        #Invert molecular problem:
        _, self.inversion_info = self.inverter.invert(ntarget, vs0, phi0, e0, ispin)


        if self.optPartition["ENS_SPIN_SYM"] is True:
            print("Warning(vpKinetic): Verify shapes of solver/vs/us")

        #Calculate the functional derivative 
        #Of the molecualr kinetic energy via euler equation

        self.V.vt = self.inverter.get_vt()

    elif self.optPartition["kinetic_part_type"] == "none":

        #Neglect the non-additive kinetic energy
        self.V.vt = np.zeros_like(self.nf)

        for i_KS in [self.KSa, self.KSb]:
            i_KS.V.vt = np.zeros_like(i_KS.n)

    elif self.optPartition["kinetic_part_type"] == "twoorbital":
        
        n1 = self.na_frac
        n2 = self.nb_frac
        eT = 0.5 * self.grid.elap

        #Evaluate kinetic energy for integer occupation

        for i_KS in [self.KSa, self.KSb]:

            i_KS.V.vt = (eT @ i_KS.n ** 0.5) / (2 * self.grid.w @ i_KS.n**0.5)

        C = self.grid.integrate( (n1**0.5 - n2**0.5)**2 )
        phi1 = (n1**0.5 - n2**0.5) / C**0.5  
        phi2 = ((n1 + n2)/2 - phi1**2)**0.5

        raise ValueError("Two orbital method not available")

    elif self.optPartition["kinetic_part_type"] == "fixed":
        pass

    else:
        raise ValueError("Kinetic energy functional not recognized")


    if self.optPartition["kinetic_part_type"] != "fixed":
        #Calculate the vp contribution
        for i_KS in [self.KSa, self.KSb]:
            i_KS.V.vp_kin = (self.V.vt - i_KS.V.vt)
        #Remove nans
            i_KS.V.vp_kin[np.isnan(i_KS.V.vp_kin)] = 0.5 

    elif self.optPartition["kinetic_part_type"] == "twoorbital":
        raise ValueError("Two orbital method not yet implemented")

