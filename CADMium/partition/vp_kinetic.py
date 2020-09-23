"""
vp_kinetic.py
"""

import sys

from .initialguessinvert import initialguessinvert
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

        for i_spin in range(nspin):
            ntarget = self.nf[:, ispin]
            if self.pol == 2:
                ntarget = 2 * ntarget

        
        phi0, e0, vs0 = self.initialguessinvert(ispin)


        




