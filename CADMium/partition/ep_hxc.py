"""
ep_hxc.py
"""
import numpy as np

def ep_hxc(self):
    """
    calculate dft components of Ep
    """
    
    if self.optPartition.interaction_type == "dft":

        ###Hartree
        ehm = 0.5 * self.hartree.v_hartree( np.sum(self.nf, axis=1) )
        Ehm = self.grid.integrate( np.sum(self.nf, axis=1) * ehm[:,0]  )

        #Fragment hartree energies
        Ehf = self.E.Ehf # Start with zero and sum fragment Eh
        ehf = np.zeros_like( ehm )

        if self.optPartition.ab_sym is True:
            # Ensemble?
            if not self.ens:
                iks = [self.KSa]
            else:
                iks = [self.KSa, self.KSA]

            for i_KS in iks:
                if hasattr(i_KS.V, "frozen") is False or i_KS.V.frozen is False:
                    ehf += np.sum( i_KS.Q, axis=1 )[:,None] * i_KS.V.eh[:, None]
                ehf += self.grid.mirror(ehf)

        else:
            if not self.ens:
                iks = [self.KSa, self.KSb]
            else:
                iks = [self.KSa, self.KSA, self.KSb, self.KSB]
            for i_KS in iks:
                if hasattr(i_KS.V, "frozen") is False or i_KS.V.frozen is False:
                    ehf += np.sum(i_KS.Q, axis=1)[:,None] * i_KS.V.eh[:, None]

        
        ###Exchange
        #epsilon_x is being brought by vp_hxc
        #Exchange energy
        Exm = self.e_x
        exm = self.epsilon_x

        #Fragment exchange energies
        Exf = self.E.Exf # Start with zero and sum fragment Eh
        exf = np.zeros_like(self.nf)

        if self.optPartition.ab_sym is True:
            if not self.ens:
                iks = [self.KSa]
            else:
                iks = [self.KSa, self.KSA]
            for i_KS in iks:
                if hasattr(i_KS.V, 'frozen') is False or i_KS.V.frozen is False:
                    exf += np.sum(i_KS.Q, axis=1)[:,None] * i_KS.V.ex
            exf += self.grid.mirror(exf)

        else:
            if not self.ens:
                iks = [self.KSa, self.KSb]
            else:
                iks = [self.KSa, self.KSA, self.KSb, self.KSB]

            for i_KS in iks:
                if hasattr(i_KS.V, 'frozen') is False or i_KS.V.frozen is False:
                    exf += np.sum(i_KS.Q, axis=1)[:,None] * i_KS.V.ex

        ###Correlation
        Ecm = self.e_c
        ecm = self.epsilon_c

        #Fragment correlation energies
        Ecf = self.E.Ecf #Start with zero and sum fragment Eh
        ecf = np.zeros_like(self.nf)

        if self.optPartition.ab_sym is True:
            if not self.ens:
                iks = [self.KSa]
            else:
                iks = [self.KSa, self.KSA]
            for i_KS in iks:
                if hasattr(i_KS.V, 'frozen') is False or i_KS.V.frozen is False:
                    ecf += np.sum(i_KS.Q, axis=1)[:, None] * i_KS.V.ec

            ecf += self.grid.mirror(ecf)

        else:
            if not self.ens:
                iks = [self.KSa, self.KSb]
            else:
                iks = [self.KSa, self.KSA, self.KSb, self.KSB]
            for i_KS in iks:
                if hasattr(i_KS.V, 'frozen') is False or i_KS.V.frozen is False:
                    ecf += np.sum(i_KS.Q, axis=1)[:, None] * i_KS.V.ec

        self.E.Ep_h = Ehm - Ehf 
        self.E.Ep_x = Exm - Exf 
        self.E.Ep_c = Ecm - Ecf 
        self.V.ep_h = ehm - ehf 
        self.V.ep_x = exm - exf 
        self.V.ep_c = ecm - ecf

    else:
        self.E.Ep_h = 0
        self.E.Ep_x = 0 
        self.E.Ep_c = 0


    ###
    if self.optPartition.hxc_part_type == "exact":
        #Partition denstiy functional theory
        self.E.Ep_hxc = self.E.Ep_h + self.E.Ep_x + self.E.Ep_c

    elif self.optPartition.hxc_part_type == "hartree":
        #Hartee only
        self.E.Ep_hxc = self.E.Ep_h

    elif self.optPartition.hxc_part_type == "overlap_xc":
        #Overlap approximation for H2
        self.vp_overlap()
        self.E.Ep_hxc = self.E.Ep_h + self.E.F * (self.E.Ep_x + self.E.Ep_c)

    elif self.optPartition.hxc_part_type == "overlap_hxc":
        "Overlap aprpoximation for H2"
        self.vp_overlap()
        self.EnsCorHar()
        self.E.Ep_hxc = self.E.F * ( self.E.Ep_h + self.E.Ep_x + self.E.Ep_c ) \
                        + (1 - self.E.F) * self.E.Ehcor

    elif self.optPartition.hxc_part_type == "overlap_hxc_2":
        self.vp_overlap()
        self.EnsCorHar()
        self.E.Ep_hxc = self.E.F * ( self.E.Ep_h + self.E.Ep_x + self.E.Ep_c ) \
                        + (1 - self.E.F) * self.E.Ehcor

    elif self.optPartition.hxc_part_type == "surprise": 
        self.vp_surprise()
        self.E.Ep_hxc = self.E.Ep_h \
                      + self.grid.integrate( self.V.s * (self.V.ep_x + self.V.ep_c) * np.sum(self.nf, axis=1))

    else: 
        raise ValueError("hxc_part_type not recognized")





