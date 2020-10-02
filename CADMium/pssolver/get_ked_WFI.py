"""
get_ked_WFI.py
"""

import numpy as np

def get_ked_WFI(self):
    """
    Collect kinetic energy densities
    """
    ked = np.zeros((self.grid.Nelem, self.pol))
    for i in range(self.pol):
        #Lets assume shape 0 is 1
        if self.solver[0,i].ked_WFI is not None:
            ked[:,i] = np.sum( self.solver[0,:].ked_WFI, axis=1 )

    return ked


    