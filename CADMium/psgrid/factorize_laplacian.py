"""

Factorize_laplacian.py
Factorizes Laplacian for Hartree calculations

"""

from scipy.sparse.linalg import splu

def factorize_laplacian(self, DISP):

    if DISP is True:
        print(" Factorizing Laplacian ... \n")

        LU = splu(self.elap, permc_spec="NATURAL")

        self.L_lap = LU.L
        self.U_lap = LU.U

        print(self.L_lap)
        print(self.U_lap)
        

    if DISP is True:
        print(" Done")