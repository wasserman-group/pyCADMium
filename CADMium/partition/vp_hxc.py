"""
vp_hxc
"""

import numpy as np

def vp_hxc(self):
    """
    Calculate dft componets of vp
    """

    if not self.ens:
        iks = [self.KSa, self.KSb]
    else:
        iks = [self.KSa, self.KSA, self.KSb, self.KSB]
    
    if self.optPartition.hxc_part_type != "nohxc" and \
       (self.optPartition.hxc_part_type == "exact" and \
        self.optPartition.interaction_type == "ni") is False:
        #Do we need to calculate non-additive dft potentials?

        #Calculate hxc functional for promolecular density
        self.V.vh = self.hartree.v_hartree(self.nf)
        self.e_x, self.epsilon_x, self.V.vx = self.exchange.get_xc(self.nf, return_epsilon=True)
        self.e_c, self.epsilon_c, self.V.vc = self.correlation.get_xc(self.nf, return_epsilon=True)

        for i_KS in iks:
            i_KS.V.vp_h = self.V.vh - i_KS.V.vh
            i_KS.V.vp_x = self.V.vx - i_KS.V.vx
            i_KS.V.vp_c = self.V.vc - i_KS.V.vc

            i_KS.V.vp_hxc = i_KS.V.vp_h + i_KS.V.vp_x + i_KS.V.vp_c

    else:
        #If not then set non-additive dft potentials to zeros
            for i_KS in iks:
                i_KS.V.vp_h = np.zeros_like(self.nf)
                i_KS.V.vp_x = np.zeros_like(self.nf)
                i_KS.V.vp_c = np.zeros_like(self.nf)
            
    
    #Do calculations for overlap type functionals
    if self.optPartition.hxc_part_type == "overlap_xc" or \
       self.optPartition.hxc_part_type == "overlap_hxc" or \
       self.optPartition.hxc_part_type == "overlap_exx":
    
        self.vp_overlap()

    #Do calculatiosn for surprisal type functionals
    if self.optPartition.hxc_part_type == "surprisal":
        pass

    #Collect hxc terms according to settings
    if self.optPartition.hxc_part_type == "exact":
        #Exactly matching
        for i_KS in iks:
            i_KS.V.vp_hxc = i_KS.V.vp_h + i_KS.V.vp_x + i_KS.V.vp_c
        
    elif self.optPartition.hxc_part_type == "hartree":
        #Hartree only
        for i_KS in iks:
            i_KS.V.vp_hxc = i_KS.V.vp_h

    elif self.optPartition.hxc_part_type == "overlap_xc":
        #Overlap approximation for H2
        #Chain rule to evaluate overlap xc term
        self.energy()
        
        for i_KS in iks:
            i_KS.V.vp_hxc =   i_KS.V.vp_h \
                            + self.E.F * (i_KS.V.vp_x + i_KS.V.vp_c) \
                            + ( self.E.Ep_x + self.E.Ep_c ) * i_KS.V.dFdn

    elif self.optPartition.hxc_part_type == "overlap_hxc":
        #Overlap approximation for H2p
        self.EnsCorHar()
        #Chain rule to evaluate overlap hxc term 
        self.energy()

        for i_KS in iks:
            i_KS.V.vp_hxc = self.E.F * (i_KS.V.vp_h + i_KS.V.vp_x + i_KS.V.vp_c) + \
                            + (self.E.Ep_h + self.E.Ep_x + self.E.Ep_c) * i_KS.V.dFdn \
                            + (1 - self.E.F) * (i_KS.V.vhcor) \
                            + (self.E.Ehcor) * i_KS.V.dFdn

    # elif self.optPartition.hxc_part_type == "overlap_hxc_2":
    #     #Overlap approximation for H2p
    #     self.EnsCorHar()
    #     #Chain rule to evaluate overlap hxc term 
    #     self.energy()
    #     for i_KS in iks:
    #         i_KS.V.vp_hxc = self.E.F * (i_KS.V.vp_h + i_KS.V.vp_x + i_KS.V.vp_c) + \
    #                         + (self.E.Ep_h + self.E.Ep_x + self.E.Ep_c) * i_KS.V.dFdn \
    #                         + (1 - self.E.F) * (i_KS.V.vhcor) \
    #                         + (self.E.Ehcor) * i_KS.V.dFdn

    elif self.optPartition.hxc_part_type == "surprisal":
        #Hmmm (?)
        self.Ep_hxc
        self.vp_surprise

        for i_KS in iks:
            i_KS.V.vp_hxc = i_KS.V.vp_h + self.V.s * np.ones((1,2)) * (i_KS.V.vp_x + i_KS.V.vp_c) \
                            + (self.V.ep_x * np.ones((1,2))) \
                            + self.V.ep_c * np.ones((1,2)) * self.nf * i_KS.V.v_s

    else:
        raise ValueError (f"{self.optPartition.hxc_part_type} is not an available hxc part method") 

                        