"""
square.py

"""
import numpy as np

def square(self, fin):
    Nspin = len(fin) / self.Nelem
    return np.reshape(fin, (self.Na, self.Nr, Nspin))