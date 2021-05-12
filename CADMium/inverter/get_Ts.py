"""
get_Ts.py
"""

import numpy as np

def get_Ts(self):
    """
    get total kinetic energy
    """

    Ts = 0.0

    if self.optInv.ens_spin_sym is not True:
        for i in range(self.solver.shape[0]):
            for j in range(self.solver.shape[1]):
                self.solver[i,j].calc_energy()
                Ts += np.sum( self.solver[i,j].Ts )

    else:
        for i in range(self.solver.shape[0]):
            self.solver[i,0].calc_energy()
            Ts += 2 * np.sum( self.solver[i,0].Ts )

    return Ts
                