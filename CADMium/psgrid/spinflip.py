"""
spinflip.py
"""

import warnings

def spinflip(self, fin):
    """
    flip spins
    """

    if len(fin.shape) == 1:
        warnings.warn('Cannot flip spin unpolarized spin density')
        return fin 
    else:
        return np.flip(fin, axis=1)