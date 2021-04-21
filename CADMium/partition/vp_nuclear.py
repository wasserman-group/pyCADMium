"""
vp_nuclear.py
"""

def vp_nuclear(self):
    """
    Calculate nuclear potential energy components of vp
    """
    
    self.KSa.V.vp_pot = self.V.vext - self.KSa.vnuc
    self.KSb.V.vp_pot = self.V.vext - self.KSb.vnuc

    if self.ens:
        self.KSA.V.vp_pot = self.V.vext - self.KSA.vnuc
        self.KSB.V.vp_pot = self.V.vext - self.KSB.vnuc