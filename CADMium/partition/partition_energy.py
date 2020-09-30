"""
partition_energy.py
""" 

def partition_energy(self):
    """
    Calculates the partiton energy
    """

    if self.optPartition["ISOLATED"] is not True:

        #Calculate components of the partition energy
        self.Ep_nuclear() #Potential energy
        self.Ep_kinetic() #Kinetic energy
        self.Ep_hxc()     #Hxc energy

        self.E.Ep = self.E.Ep_pot + self.E.Ep_kin + self.E.Ep_hxc


    else:
        self.E.Ep     = 0.0
        self.E.Ep_pot = 0.0
        self.E.Ep_kin = 0.0
        self.E.Ep_hxc = 0.0
    