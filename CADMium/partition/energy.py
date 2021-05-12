"""
energy.py
"""

from copy import copy as copy

def energy(self):
    """
    Calculate energies
    """

    #Sum of fragment energies
    Tsa       = 0.0
    Eksa      = 0.0
    Enuca     = 0.0
    Exa       = 0.0
    Eca       = 0.0
    Eha       = 0.0
    Vhxca     = 0.0
    self.E.Ea = 0.0

    for i_KS in [self.KSa]:
        if i_KS.E is None:
            i_KS.energy()

        self.E.Ea += i_KS.scale * i_KS.E.E
        Tsa       += i_KS.scale * i_KS.E.Ts
        Eksa      += i_KS.scale * i_KS.E.Eks
        Enuca     += i_KS.scale * i_KS.E.Enuc

        if self.optPartition.interaction_type == "dft":
            Exa   += i_KS.scale * i_KS.E.Ex
            Eca   += i_KS.scale * i_KS.E.Ec
            Eha   += i_KS.scale * i_KS.E.Eh
            Vhxca += i_KS.scale * i_KS.E.Vhxc

    if not self.optPartition.ab_sym:
        Tsb       = 0.0
        Eksb      = 0.0
        Enucb     = 0.0
        Exb       = 0.0
        Ecb       = 0.0
        Ehb       = 0.0
        Vhxcb     = 0.0
        self.E.Eb = 0.0

        for i_KS in [self.KSb]:
            if i_KS.E is None:
                i_KS.energy()

            self.E.Eb += i_KS.scale * i_KS.E.E
            Tsb       += i_KS.scale * i_KS.E.Ts
            Eksb      += i_KS.scale * i_KS.E.Eks
            Enucb     += i_KS.scale * i_KS.E.Enuc

            if self.optPartition.interaction_type == "dft":
                Exb   += i_KS.scale * i_KS.E.Ex
                Ecb   += i_KS.scale * i_KS.E.Ec
                Ehb   += i_KS.scale * i_KS.E.Eh
                Vhxcb += i_KS.scale * i_KS.E.Vhxc

    else:

        Tsb       = copy(Tsa)
        Eksb      = copy(Eksa)
        Enucb     = copy(Enuca)
        Exb       = copy(Exa)
        Ecb       = copy(Eca)
        Ehb       = copy(Eha)
        Vhxcb     = copy(Vhxca)
        self.E.Eb = copy(self.E.Ea)    

    if self.optPartition.ens_spin_sym:
        #We double the fragment energies to account
        #for the spin flipped components

        Tsa       *= 2.0
        Eksa      *= 2.0
        Enuca     *= 2.0
        Exa       *= 2.0
        Eca       *= 2.0
        Eha       *= 2.0
        Vhxca     *= 2.0
        self.E.Ea *= 2.0

        Tsb       *= 2.0
        Eksb      *= 2.0
        Enucb     *= 2.0
        Exb       *= 2.0
        Ecb       *= 2.0
        Ehb       *= 2.0
        Vhxcb     *= 2.0
        self.E.Eb *= 2.0

    #Sum of fragment energies
    self.E.Ef    = self.E.Ea + self.E.Eb
    self.E.Tsf   = Tsa + Tsb
    self.E.Eksf  = Eksa + Eksb
    self.E.Enucf = Enuca + Enucb
    self.E.Exf   = Exa + Exb
    self.E.Ecf   = Eca + Ecb
    self.E.Ehf   = Eha + Ehb
    self.E.Vhxcf = Vhxca + Vhxcb

    #Calculate partition energy
    self.partition_energy()

    #Total electronic energy
    self.E.Et = self.E.Ef + self.E.Ep

    #Nuclear nuclear repulsion
    self.E.Vnn = self.Za * self.Zb / (2 * self.grid.a)

    self.E.E = self.E.Et + self.E.Vnn

    self.E.evals_a = []
    self.E.evals_b = []

    if self.optPartition.ab_sym is True:
        for i_KS in [self.KSa]:
            self.E.evals_a = i_KS.E.evals
        self.E.evals_b = self.E.evals_a
    
    else:
        for i_KS in [self.KSa, self.KSb]:
            self.E.evals_a = i_KS.E.evals
            self.E.evals_b = i_KS.E.evals



    
    

