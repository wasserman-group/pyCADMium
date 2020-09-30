"""
ep_kinetic.py
"""

import numpy as np
import sys

def ep_kinetic(self):
    """
    Calculate kinetic energy components of Ep
    """

    #Select method for calculating kinetic energy

    if self.optPartition["kinetic_part_type"] == "vonweiz":

        #Use the von-Weizacker inversion
        Tsm =   (2 * np.pi * self.grid.hr * self.grid.ha ) \
              * np.sum(self.grid.wi * np.sum(self.nf, axis=1)**0.5  \ 
              * (0.5 * self.grid.elap @ np.sum(self.nf, axis=1)**0.5))

        tsm =   (np.sum(self.nf, axis=1)**(-0.5)) / self.grid.w \
              * (-0.5 * sefl.grid.elap @ np.sum(self.nf, axis=1)**0.5 )


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


    elif self.optPartition["kinetic_part_type"] == "libxcke" or \
         self.optPartition["kinetic_part_{ype"] == "paramke":
         #Use a kinetic energy from libxc

        raise ValueError(f"Kinetic_part_type {self.optPartition['kinetic_part_type']} has not been implemented")


    elif self.optPartition["inversion"] == "inversion":
        pass

    print("leaving from ep_kinetic")
    sys.exit()