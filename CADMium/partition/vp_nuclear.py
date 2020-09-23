"""
vp_nuclear.py
"""

def vp_nuclear(self):
    """
    Calculate nuclear potential energy components of vp
    """

    self.KSa.V.vpot = self.V.vext - self.KSa.vnuc
    self.KSb.V.vpot = self.V.vext - self.KSb.vnuc
