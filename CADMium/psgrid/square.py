"""
square.py

"""
import numpy as np

def square(self, fin):
    Nspin = int(len(fin) / self.Nelem)
    square = np.reshape(fin, (self.Na, self.Nr, Nspin)) 
    return square