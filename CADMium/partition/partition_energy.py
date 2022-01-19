"""
partition_energy.py
""" 

def partition_energy(self):
    """
    Calculates the partiton energy
    """

    if not self.optPartition.isolated:
        # Calculate components of the partition energy
        self.ep_nuclear() #Potential energy
        self.ep_kinetic() #Kinetic energy
        self.ep_hxc()     #Hxc energy
        self.E.Ep = self.E.Ep_pot + self.E.Ep_kin + self.E.Ep_hxc

    else:
        self.E.Ep     = 0.0
        self.E.Ep_pot = 0.0
        self.E.Ep_kin = 0.0
        self.E.Ep_hxc = 0.0
    