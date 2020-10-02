"""
get_Ts.py
"""

import numpy as np

def get_Ts(self):
    """
    get total kinetic energy
    """

    Ts = 0.0

    if self.optInversion["ENS_SPIN_SYM"] is not True:
        
        for i in range(self.solver.shape[0]):
            for j in range(self.solver.shape[1]):
                #Get_ts cannot be computed as a pssolver method
                #Ts += np.sum( self.solver[i,j].get_Ts() )
                self.solver[i,j].calc_energy()
                Ts += np.sum( self.solver[i,j].Ts )

    else:

        for i in range(self.solver.shape[0]):
            self.solver[i,0].calc_energy()
            Ts += 2 * np.sum( self.solver[i,0].Ts )

    return Ts
                