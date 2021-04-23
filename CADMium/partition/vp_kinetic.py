"""
vp_kinetic.py
"""

import numpy as np
# np.set_printoptions(precision=8)

def vp_kinetic(self):
    """
    Calculate Kinetic energy components of vp
    """

    #---> Select method for calculating vp Kinetic
    if self.optPartition.kinetic_part_type == "vonweiz":

        # Use the von-weizacker inversion
        self.V.vt = (-0.5 * self.grid.elap @ (self.nf**0.5)) / (self.nf**0.5) / (self.grid.w * np.ones_like(self.nf).T).T

        if not self.ens:
            iksa = [self.KSa]
            iksb = [self.KSb]
        else:
            iksa = [self.KSa, self.KSA]
            iksb = [self.KSb, self.KSB]

        for KSa in iksa: 
            # Evaluate kinetic energy for integer ocupations
            KSa.V.vt = (-0.5 * self.grid.elap @ (KSa.n**0.5)) / (KSa.n**0.5) / np.repeat(self.grid.w[:,None], KSa.n.shape[1], axis=1)
            # Get rid of infties/nans  
            KSa.V.vt = np.nan_to_num(KSa.V.vt, nan=0.0, posinf=0.0, neginf=0.0)

        for KSb in iksb:   
            KSb.V.vt = (-0.5 * self.grid.elap @ (KSb.n**0.5)) / (KSb.n**0.5) / np.repeat(self.grid.w[:,None], KSb.n.shape[1], axis=1)
            KSb.V.vt = np.nan_to_num(KSb.V.vt, nan=0.0, posinf=0.0, neginf=0.0)

    elif self.optPartition.kinetic_part_type == "libxcke" or \
         self.optPartition.kinetic_part_type == "paramke":
        
        #Use a kinetic energy functional from libxc
        #Evaluate MOLECULAR kinetic energy
        self.Tsm, self.tsm, self.V.vt = self.kinetic.get_xc(self.nf, return_epsilon=True)

        if not self.ens:
            iks = [self.KSa, self.KSb]
        else:
            iks = [self.KSa, self.KSA, self.KSb, self.KSB]

        # Evaluate FRAGMENT kinetic energy
        for iKS in iks:
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

        if self.ens:
            for i in range(self.Nmo_A.shape[0]):
                for j in range(self.Nmo_A.shape[1]):
                    if self.KSA.solver[i,j].homo is not None:
                        ua = np.max( (ua, self.KSA.solver[i,j].homo) )
                    if self.KSB.solver[i,j].homo is not None:
                        ua = np.max( (ub, self.KSB.solver[i,j].homo) )
        u = np.max((ua, ub))
        nspin = self.pol - 1

        if self.optPartition.ens_spin_sym:
            u = max([u])
            nspin = 0

        self.KSa.V.vt = np.full(self.KSa.veff.shape, u) - self.KSa.veff
        self.KSb.V.vt = np.full(self.KSb.veff.shape, u) - self.KSb.veff
        if self.ens:
            self.KSA.V.vt = np.full(self.KSA.veff.shape, u) - self.KSA.veff
            self.KSB.V.vt = np.full(self.KSB.veff.shape, u) - self.KSB.veff

        for ispin in range(nspin+1):
            ntarget = self.nf[:, ispin]
            if self.pol == 2:
                ntarget = 2 * ntarget
            phi0, e0, vs0 = self.initialguessinvert(ispin)
            phi0 = phi0.flatten()[:,None]
            e0   = e0.flatten()

            #Invert molecular problem:
            _, self.inversion_info = self.inverter.invert(ntarget, vs0, phi0, e0, ispin)

        # print("Leaving through vp_kinetic")
        # sys.exit()

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
        if not self.ens:
            iks = [self.KSa, self.KSb]
        else:
            iks = [self.KSa, self.KSA, self.KSb, self.KSB]
        for i_KS in iks:
            i_KS.V.vt = np.zeros_like(i_KS.n)

    elif self.optPartition.kinetic_part_type == "twoorbital":
        
        n1 = self.na_frac
        n2 = self.nb_frac
        eT = 0.5 * self.grid.elap

        #Evaluate kinetic energy for integer occupation
        if not self.ens:
            iks = [self.KSa, self.KSb]
        else:
            iks = [self.KSa, self.KSA, self.KSb, self.KSB]

        for i_KS in iks:
            i_KS.V.vt = (eT @ i_KS.n ** 0.5) / (2 * self.grid.w @ i_KS.n**0.5)

        C = self.grid.integrate( (n1**0.5 - n2**0.5)**2 )
        phi1 = (n1**0.5 - n2**0.5) / C**0.5  
        phi2 = ((n1 + n2)/2 - phi1**2)**0.5

        raise ValueError("Two orbital method completed yet")

    elif self.optPartition.kinetic_part_type == "fixed":
        pass

    else:
        raise ValueError("Kinetic energy functional not recognized")

    #---> Finalize VP Kinetic
    if self.optPartition.kinetic_part_type == "twoorbital":
        """
        Two orbital requires special treatment because
        molecular funcitonal derivative depens on the fragment
        """

        i = 0
        if not self.ens:
            iks = [self.KSa, self.KSb]
        else:
            iks = [self.KSa, self.KSA, self.KSb, self.KSB]

        for IKS in iks:
            IKS.V.vp_kin = self.V.vt[:,i] - IKS.V.vt
            i += 1 

    elif self.optPartition.kinetic_part_type != "fixed":

        if not self.ens:
            iks = [self.KSa, self.KSb]
        else:
            iks = [self.KSa, self.KSA, self.KSb, self.KSB]

        #Calculate the vp contribution
        for i_KS in iks:
            i_KS.V.vp_kin = (self.V.vt - i_KS.V.vt)
        #Remove nans
            i_KS.V.vp_kin[np.isnan(i_KS.V.vp_kin)] = 0.0


