"""

mirror.py
Mirror function accros AB plane

"""
import numpy as np

def mirror(self, fin):
    fout = self.square(fin)
    return np.reshape(fout[::-1,:,:],fin.shape)
    