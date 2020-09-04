"""
reduced_grad.py
"""

import numpy as np

def reduced_grad(self, n):
    """
    Calculate the reduced density gradient
    """

    pol = n.shape[1]

    if pol == 1:
        s = self.sigma(n) ** (0.5) / (2 * (3*np.pi**2)**(1/3) @ n** (4/3) )

    if pol == 2:
        s = np.zeros_like(n)
        s[:, 0] = self.sigma(2 * n[:,0])**(0.5) / (2 * (3*np.pi**2)**(1/3) @ (2 * n[:,0])**(4/3))
        s[:, 0] = self.sigma(2 * n[:,1])**(0.5) / (2 * (3*np.pi**2)**(1/3) @ (2 * n[:,1])**(4/3))

    
    return s
