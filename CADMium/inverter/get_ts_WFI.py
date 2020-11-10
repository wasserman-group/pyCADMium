"""
get_ts_WFI.py
"""

import numpy as np

def get_ts_WFI(self):
    """
    Get kinetic energy density
    """

    

    ts = np.zeros((self.grid.Nelem, len(self.solver[0,:])  ))

    if self.optInversion["ENS_SPIN_SYM"] is not True:

        for i in range(self.solver.shape[0]):
            for j in range(self.solver.shape[1]):
                self.solver[i,j].calc_ked_WFI()
                #ts[i,j] = self.solver[i,j].get_ked_WFI()
                
        #get_ked_WFI cannot be defined as a solver's method
        #Get Kinetic Energy Density
        for i in range(self.solver.shape[0]):
            for j in range(self.solver.shape[1]):
                if self.solver[i,j].ked_WFI is not None:
                    ts[:,j] = np.sum( self.solver[i,j].ked_WFI, axis=1 ) 

    else:
        for i in range(self.solver.shape[0]):
            self.solver[i,0].calc_ked_WFI()

        #Get Kinetic Energy Density
        for i in range(self.solver.shape[0]):
            for j in range(self.solver.shape[1]):
                if self.solver[i,j].ked_WFI is not None:
                    ts[:,j] = np.sum( self.solver[i,j].ked_WFI, axis=1 ) 

    return ts