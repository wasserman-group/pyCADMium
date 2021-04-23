"""
partition_potential
"""

import numpy as np

def partition_potential(self):
    """
    Calculate the partition potential
    """
 
    if self.optPartition.vp_calc_type== "component":
        
        #Calculate components
        #Potential energy compoenent
        self.vp_nuclear()
        #Calculate kinetic term
        self.vp_kinetic()
        #Calculate hxc terms
        self.vp_hxc()


        #Build the partition potential and components
        self.V.vp     = np.zeros_like(self.nf)
        self.V.vp_pot = np.zeros_like(self.nf)
        self.V.vp_kin = np.zeros_like(self.nf)
        self.V.vp_hxc = np.zeros_like(self.nf)
        self.V.vp_h   = np.zeros_like(self.nf)
        self.V.vp_x   = np.zeros_like(self.nf)
        self.V.vp_c   = np.zeros_like(self.nf)

        if not self.ens:
            iks = [self.KSa, self.KSb]
        else:
            iks = [self.KSa, self.KSA, self.KSb, self.KSB]

        for i_KS in iks:
            i_KS.V.vp = i_KS.V.vp_pot + i_KS.V.vp_kin + i_KS.V.vp_hxc

            self.V.vp     +=  i_KS.V.vp     * i_KS.Q
            self.V.vp_pot +=  i_KS.V.vp_pot * i_KS.Q
            self.V.vp_kin +=  i_KS.V.vp_kin * i_KS.Q
            self.V.vp_hxc +=  i_KS.V.vp_hxc * i_KS.Q
            self.V.vp_h   +=  i_KS.V.vp_h   * i_KS.Q
            self.V.vp_x   +=  i_KS.V.vp_x   * i_KS.Q
            self.V.vp_c   +=  i_KS.V.vp_c   * i_KS.Q

        if self.optPartition.ens_spin_sym is True:
            self.V.vp     += self.grid.spinflip(self.V.vp)
            self.V.vp_pot += self.grid.spinflip(self.V.vp_pot)
            self.V.vp_kin += self.grid.spinflip(self.V.vp_kin)
            self.V.vp_hxc += self.grid.spinflip(self.V.vp_hxc)
            self.V.vp_h   += self.grid.spinflip(self.V.vp_h)
            self.V.vp_x   += self.grid.spinflip(self.V.vp_x)
            self.V.vp_c   += self.grid.spinflip(self.V.vp_c)

        if self.optPartition.ab_sym is True:
            self.V.vp     = 0.5 * ( self.V.vp + self.grid.mirror(self.V.vp ))
            self.V.vp_pot = 0.5 * ( self.V.vp_pot + self.grid.mirror(self.V.vp_pot) )
            self.V.vp_kin = 0.5 * ( self.V.vp_kin + self.grid.mirror(self.V.vp_kin) )
            self.V.vp_hxc = 0.5 * ( self.V.vp_hxc + self.grid.mirror(self.V.vp_hxc) )
            self.V.vp_h   = 0.5 * ( self.V.vp_h + self.grid.mirror(self.V.vp_h) )
            self.V.vp_x   = 0.5 * ( self.V.vp_x + self.grid.mirror(self.V.vp_x) )
            self.V.vp_c   = 0.5 * ( self.V.vp_c + self.grid.mirror(self.V.vp_c) )

        vp = self.V.vp

            
    elif self.optPartition.vp_calc_type == 'potential_inversion':

        #Find kinetic energy functional derivative for
        #fragments using euler equation

        u = max((self.KSa.u, self.KSb.u))

        if not self.ens:
            iks = [self.KSa, self.KSb]
        else:
            iks = [self.KSa, self.KSA, self.KSb, self.KSB]

        for i_KS in iks:
            i_KS.V.vt = u - i_KS.veff

        #Form initial guess for molecular inversion
        if self.pol == 2:
            vs0 = np.zeros_like(self.nf)
            for i_KS in iks:
                vs0 += i_KS.veff * i_KS.Q
                if hasattr(i_KS.V, 'vp_kin') is True:
                    vs0 -= i_KS.V.vp_kin * i_KS.Q
            if self.optPartition.ens_spin_sym:
                vs0 += self.grid.spinfip(vs0)
        else:
            phi0, e0, vs0 = self.initialguessinvert()

        #Invert molecular problem
        _, self.inversion_info = self.inverter.invert(self.nf, vs0, phi0, e0)

        #Calculate the functional derivative of the
        #molecular kinetic energy via the euler equation
        self.V.vt = self.inverter.get_vt()
        self.V.vs = self.inverter.vs

        #Calculate hxc functional for promolecular density
        self.V.vh = self.hartree.v_hartree(self.nf)
        _, self.V.vx = self.exchange.get_xc(self.nf)
        _, self.V.vc = self.correlation.get_xc(self.nf)
        self.V.vhxc = self.V.vh + self.V.vx + self.V.vc

        if hasattr(self.V, 'vp') is False:
            self.V.vp = np.zeros_like(self.nf)

        self.V.vp += self.Beta * (self.V.vext + self.V.vhxc - self.V.vs)

        vp = self.V.vp

    return vp
    
