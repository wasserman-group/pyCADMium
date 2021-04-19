"""
vp_kinetic.py
"""

import numpy as np
# np.set_printoptions(precision=8)

def vp_kinetic(self):
    """
    Calculate Kinetic energy components of vp
    """

    # Select method for calculating kinetic energy
    if self.optPartition.kinetic_part_type == "vonweiz":

        # Use the von-weizacker inversion
        self.V.vt = (-0.5 * self.grid.elap @ (self.nf**0.5)) / (self.nf**0.5) / (self.grid.w * np.ones_like(self.nf).T).T

        # Evaluate kinetic energy for integer ocupations
        self.KSa.V.vt = (-0.5 * self.grid.elap @ (self.KSa.n**0.5)) / (self.KSa.n**0.5) / np.repeat(self.grid.w[:,None], self.KSa.n.shape[1], axis=1)
        self.KSb.V.vt = (-0.5 * self.grid.elap @ (self.KSb.n**0.5)) / (self.KSb.n**0.5) / np.repeat(self.grid.w[:,None], self.KSb.n.shape[1], axis=1)

        # Get rid of infties/nans  
        self.KSa.V.vt = np.nan_to_num(self.KSa.V.vt, nan=0.0, posinf=0.0, neginf=0.0)
        self.KSb.V.vt = np.nan_to_num(self.KSb.V.vt, nan=0.0, posinf=0.0, neginf=0.0)

    elif self.optPartition.kinetic_part_type == "libxcke" or \
         self.optPartition.kinetic_part_type == "paramke":
        
        #Use a kinetic energy functional from libxc
        #Evaluate MOLECULAR kinetic energy
        self.Tsm, self.tsm, self.V.vt = self.kinetic.get_xc(self.nf, return_epsilon=True)

        # Evaluate FRAGMENT kinetic energy
        for iKS in [self.KSa, self.KSb]:
            iKS.V.Ts, iKS.V.ts, iKS.V.vt = self.kinetic.get_xc(iKS.n, return_epsilon=True)            

    elif self.optPartition.kinetic_part_type == "inversion":
        #Find kinetic energy fucntional derivative for fragmetns using Euler equation
        #Fragments using euler equations
        # u = max(self.KSa.u, self.KSb.u)
        ua, ub = -1 / np.spacing(1), -1 / np.spacing(1)

        #Check the values of each solver's homos and pick the largest
        for i in range(self.Nmo_a.shape[0]):
            for j in range(self.Nmo_a.shape[1]):
                if self.KSa.solver[i,j].homo is not None:
                    ua = np.max( (ua, self.KSa.solver[i,j].homo) )
                if self.KSb.solver[i,j].homo is not None:
                    ua = np.max( (ub, self.KSb.solver[i,j].homo) )
        u = np.max((ua, ub))
        nspin = self.pol - 1

        if self.optPartition.ens_spin_sym:
            u = max([u])
            nspin = 0

        self.KSa.V.vt = np.full(self.KSa.veff.shape, u) - self.KSa.veff
        self.KSb.V.vt = np.full(self.KSb.veff.shape, u) - self.KSb.veff

        for ispin in range(nspin+1):
            ntarget = self.nf[:, ispin]
            if self.pol == 2:
                ntarget = 2 * ntarget

            phi0, e0, vs0 = self.initialguessinvert(ispin)

            #Invert molecular problem:
            _, self.inversion_info = self.inverter.invert(ntarget, vs0, phi0, e0, ispin)

        if self.optPartition.ens_spin_sym is True:
            self.inverter.solver[:,1] = self.inverter.solver[:,0]
            self.inverter.vs[:,1] = self.inverter.vs[:,0]
            self.inverter.us[:,1] = self.inverter.us[:,0]

        #Calculate the functional derivative 
        #Of the molecualr kinetic energy via euler equation
        self.V.vt = self.inverter.get_vt()

    elif self.optPartition.kinetic_part_type == "none":

        #Neglect the non-additive kinetic energy
        self.V.vt = np.zeros_like(self.nf)

        for i_KS in [self.KSa, self.KSb]:
            i_KS.V.vt = np.zeros_like(i_KS.n)

    elif self.optPartition.kinetic_part_type == "twoorbital":
        
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

    elif self.optPartition.kinetic_part_type == "fixed":
        pass

    else:
        raise ValueError("Kinetic energy functional not recognized")

    if self.optPartition.kinetic_part_type != "fixed":
        #Calculate the vp contribution
        for i_KS in [self.KSa, self.KSb]:
            i_KS.V.vp_kin = (self.V.vt - i_KS.V.vt)
        #Remove nans
            i_KS.V.vp_kin[np.isnan(i_KS.V.vp_kin)] = 0.0

    elif self.optPartition.kinetic_part_type == "twoorbital":
        raise ValueError("Two orbital method not yet implemented")

