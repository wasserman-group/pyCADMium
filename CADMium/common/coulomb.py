"""
coulomb.py
"""

import numpy as np

def coulomb(grid, Za, Zb):
    """
    Calculates the external field from nuclei with charges Za, Zb
    """

    v = -1.0/grid.a * ((Za + Zb) * np.cosh(grid.Xr) - (Za - Zb) * np.cos(grid.Xa))  / (
                                       np.cosh(grid.Xr)**2 - np.cos(grid.Xa)**2
                                       )
    
    return v
                                        