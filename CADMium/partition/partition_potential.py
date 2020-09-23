"""
partition_potential
"""

import numpy as np

def partition_potential(self):
    """
    Calculate the partition potential
    """
 
    if self.optPartition["vp_calc_type"] == "component":

        #Calculate components
        self.vp_nuclear()
        self.vp_kinetic()
        # self.vp_hxc()

        
        #Build the partition potential and components

        vp = np.zeros([1])
    return vp