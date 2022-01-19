"""
ep_kinetic.py
"""
import sys
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

        if not self.ens:
            iks = [self.KSa, self.KSb]
        else:
            iks = [self.KSa, self.KSA, self.KSb, self.KSB]

        for i_KS in iks:
            Tsf += i_KS.scale * ( 2 * np.pi * self.grid.hr * self.grid.ha ) \
                              * np.sum( self.grid.wi * np.sum(i_KS.n, axis=1)**0.5 ) \
                              * ( -0.5 * self.grid.elap @ np.sum(i_KS.n, axis=1)**0.5 )

            i_KS.V.ts = np.sum( i_KS.n, axis=1 )**(-0.5) / self.grid.w \
                        * (-0.5 * self.grid.elap @ np.sum( i_KS.n, axis=1 )**0.5)

            tsf += i_KS.scale * np.sum(i_KS.Q, axis=1) * i_KS.V.ts

    elif self.optPartition.kinetic_part_type == "libxcke" or \
         self.optPartition.kinetic_part_type == "paramke":
         #Use a kinetic energy from libxc

        # Evaluate MOLECULAR kinetic energy
        # Tsm and tsm calculated in ep_kinetic
        tsm = self.tsm[:,0]
        Tsm = self.Tsm

        #Evaluate kinetic energy for integer occupation
        Tsf = 0 #Start with zero and sum fragment Kinetic energy
        tsf = np.zeros_like( tsm )

        if not self.ens:
            iks = [self.KSa, self.KSb]
        else:
            iks = [self.KSa, self.KSA, self.KSb, self.KSB]

        # Ts and ts for fragments obtained in ep_kinetic
        for iKS in iks:
            Tsf += iKS.scale * iKS.V.Ts
            tsf += iKS.scale * np.sum( iKS.Q, axis=1 ) * iKS.V.ts[:,0]

    elif self.optPartition.kinetic_part_type == "inversion":
        
        tsm = np.sum( self.inverter.get_ts_WFI(), axis=1 ) / np.sum( self.nf, axis=1 )
        Tsm = self.inverter.get_Ts()

        #Evaluate kinetic energy for fragments
        Tsf = 0.0
        tsf = np.zeros_like(tsm)
        tsf = tsf[:, None]

        if self.optPartition.ab_sym is not True:
            if not self.ens:
                KSab = [self.KSa, self.KSb]
            else:
                KSab = [self.KSa, self.KSA, self.KSb, self.KSB]
        else:
            if not self.ens:
                KSab = [self.KSa]
            else:
                KSab = [self.KSa, self.KSA]

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

        n1 = self.na_frac
        n2 = self.nb_frac
        eT = -1/2 * self.grid.elap 

        integrand = (n1**(0.5) - n2**(0.5))**2
        C  = self.grid.integrate( np.sum(integrand, axis=1) )

        phi1 = ( (self.grid.integrate( np.sum(n1+n2, axis=1) ))/2 )**(0.5) * (n1**(0.5) - n2**(0.5))/C**(0.5)  
        phi2 = ((n1+n2) - phi1**2)**(0.5) 

        TO_phi1 = self.grid.integrate( (phi1 * (eT @ phi1) / self.grid.w[:,None])[:,0] )
        TO_phi2 = self.grid.integrate( (phi2 * (eT @ phi2) / self.grid.w[:,None])[:,0] )

        self.V.phi1 = phi1
        self.V.phi2 = phi2

        Tsm = TO_phi1 + TO_phi2
        tsm = (phi1 * (eT @ phi1) + phi2 * (eT @ phi2)) / self.grid.w[:,None] / np.sum(self.nf, axis=1)[:,None]
        self.V.tsm = tsm

        # Evaluate kinetic energy for fragments
        Tsf = 0 # Start with zero and sum fragment kinetic energy. 
        tsf = np.zeros_like( tsm )

        if not self.optPartition.ab_sym:
            KSab = [self.KSa, self.KSb]
        else:
            KSab = [self.KSa]
        
        for i_KS in KSab:
            
            Tsf +=  i_KS.scale * ( 2 * np.pi * self.grid.hr * self.grid.ha ) * np.sum( self.grid.wi * (np.sum(i_KS.n, axis=1)**0.5) * (eT@np.sum(i_KS.n, axis=1)**0.5) ) 
            i_KS.V.ts = i_KS.scale * ( np.sum(i_KS.n, axis=1)**(0.5) ) * (eT @ np.sum( i_KS.n, axis=1 )**(0.5))  / (self.grid.w * np.sum(i_KS.n, axis=1) )

            # print("1", ( np.sum(i_KS.n, axis=1)**(0.5) ) )
            # print("2", (eT @ np.sum( i_KS.n, axis=1 )**(0.5)) )
            # print("3", (self.grid.w * np.sum(i_KS.n, axis=1) ) )

            tsf += i_KS.Q * i_KS.V.ts[:,None]

        if self.optPartition.ab_sym:
            tsf += tsf + self.grid.mirror(tsf)
            Tsf *= 2      

        # print("Yo soy Tsf", tsf)
        # sys.exit()


    elif self.optPartition.kinetic_part_type == "fixed":
        raise ValueError("Fixed method not yet implemented")

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

        if not self.ens:
            KSab = [self.KSa, self.Ksb]
        else:
            KSab = [self.KSa, self.KSA, self.KSb, self.KSB]

        Ep_kin = self.Ep_kin0

        for i_KS in KSab:
            if i_KS.frozen is not True:
                Ep_kin += self.grid.integrate( np.sum( i_KS.vp_kin * i_KS.n, axis=1 ) )

        self.E.Ep_kin = Ep_kin
    
