"""
ep_kinetic.py
"""

import numpy as np

def ep_kinetic(self):
    """
    Calculate kinetic energy components of Ep
    """

    #Select method for calculating kinetic energy

    if self.optPartition.kinetic_part_type == "vonweiz":

        #Use the von-Weizacker inversion
        Tsm =   (2 * np.pi * self.grid.hr * self.grid.ha ) \
              * np.sum(self.grid.wi * np.sum(self.nf, axis=1)**0.5 \
              * (0.5 * self.grid.elap @ np.sum(self.nf, axis=1)**0.5))

        tsm =   (np.sum(self.nf, axis=1)**(-0.5)) / self.grid.w \
              * (-0.5 * self.grid.elap @ np.sum(self.nf, axis=1)**0.5 )

        #Evaluate kinetic energy for integer ocupation 
        #Densities
        Tsf = 0.0 #Start with zero and sum fragment kinetic energy
        tsf = np.zeros_like(tsm)

        for i_KS in [self.KSa, self.KSb]:
            Tsf += i_KS.scale * ( 2 * np.pi * self.grid.hr * self.grid.ha ) \
                              * np.sum( self.grid.wi * np.sum(i_KS.n, axis=1)**0.5 ) \
                              * ( -0.5 * self.grid.elap @ np.sum(i_KS.n, axis=1)**0.5 )

            i_KS.V.ts = np.sum( i_KS.n, axis=1 )**(-0.5) / self.grid.w \
                        * (-0.5 * self.grid.elap @ np.sum( i_KS.n, axis=1 )**0.5)

            tsf += i_KS.scale * np.sum(i_KS.Q, axis=1) * i_KS.V.ts


    elif self.optPartition.kinetic_part_type == "libxcke" or \
         self.optPartition.kinetic_part_type == "paramke":
         #Use a kinetic energy from libxc

        raise ValueError(f"Kinetic_part_type {self.optPartition.kinetic_part_type} \
                           has not been implemented")

    elif self.optPartition.kinetic_part_type == "inversion":
        
        tsm = np.sum( self.inverter.get_ts_WFI(), axis=1 ) / np.sum( self.nf, axis=1 )
        Tsm = self.inverter.get_Ts()

        #Evaluate kinetic energy for fragments
        Tsf = 0.0
        tsf = np.zeros_like(tsm)
        tsf = tsf[:, None]

        if self.optPartition.ab_sym is not True:
            KSab = [self.KSa, self.KSb]

        else:
            KSab = [self.KSa]

        for i_KS in KSab:
            assert i_KS.E.Ts != None, ("Fragment energies must be evaluated before partition energy")
            Tsf += i_KS.scale * i_KS.E.Ts
            
            i_KS.V.ts = np.zeros_like( i_KS.n )
            ts = np.zeros((self.grid.Nelem, self.pol))

            for i in range(i_KS.solver.shape[0]):
                for j in range(i_KS.solver.shape[1]):
                    i_KS.solver[i,j].calc_ked_WFI()
                    
            #Same as running get_ked_WFI for solver object
            for i in range(i_KS.solver.shape[0]):
                for j in range(i_KS.solver.shape[1]):
                    if i_KS.solver[i,j].ked_WFI is not None:
                        ts[:,j] = np.sum( i_KS.solver[i,j].ked_WFI, axis=1 ) 

            for i in range(i_KS.solver.shape[0]):
                for j in range(i_KS.solver.shape[1]):
                    i_KS.V.ts += np.sum( ts, axis=1)[:,None] / np.sum(i_KS.n, axis=1)[:,None]

            tsf += i_KS.scale * np.sum( i_KS.Q, axis=1 )[:,None] * np.sum( i_KS.V.ts, axis=1 )[:,None]

        if self.optPartition.ab_sym is True:
            tsf += self.grid.mirror(tsf)
            Tsf *= 2.0

    elif self.optPartition.kinetic_part_type == "none":
        Tsm = 0.0
        tsm = np.zeros((self.grid.Nelem, 1))
        
        Tsf = 0.0
        tsf = np.zeros((self.grid.Nelem, 1))

    elif self.optPartition.kinetic_part_type == "twoorbital":
        raise ValueError("Twoorbital method not yet implemented")

    elif self.optPartition.kinetic_part_type == "fixed":
        pass

    else:
        raise ValueError("Kinetic energy method not recognized")


    ###
    if self.optPartition.kinetic_part_type != "fixed":
        #if ENS_SPIN_SYM set we double the fragment energies to
        #account for the spin flipped components

        if self.optPartition.ens_spin_sym:
            Tsf *= 2.0
            tsf *= 2.0

        self.E.Ep_kin = Tsm - Tsf
        self.V.ep_kin = tsm - tsf

    else:
        #Fixed vp_kinetic
        KSab = [self.KSa, self.Ksb]
        Ep_kin = self.Ep_kin0

        for i_KS in KSab:
            if i_KS.frozen is not True:
                Ep_kin += self.grid.integrate( np.sum( i_KS.vp_kin * i_KS.n, axis=1 ) )

        self.E.Ep_kin = Ep_kin
    
