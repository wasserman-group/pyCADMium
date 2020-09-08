"""
calc_nuclear_potential.py
"""

from ..common.coulomb import coulomb

def calc_nuclear_potential(self):

    self.va = coulomb(self.grid, self.Za, 0)
    self.vb = coulomb(self.grid, 0, self.Zb)

    self.V.vext = self.va[:, np.ones(1, self.pol)] + \
                  self.vb[:, np.ones(1, self.pol)]