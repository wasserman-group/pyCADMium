"""
normalize_orbitals.py
"""

def normalize_orbitals(self):
    
    #Figure out occupation number
    if self.m == 0:
        occ = 2
    else:
        occ = 4

    if self.pol == 2:
        occ = occ / 2

    for i in range(self.Nmo):
        self.phi[:,i] = occ**(0.5) * self.phi[:,i] / self.grid.integrate(self.phi[:,i]**2)**(0.5)