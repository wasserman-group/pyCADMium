"""
square.py

"""
import numpy as np

def square(self, fin):
    Nspin = int(fin.size / self.Nelem)
    square = np.reshape(fin, (self.Na, self.Nr, Nspin), order="F") 
    return square