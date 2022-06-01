import numpy as np

from CADMium import Psgrid
from CADMium import Kohnsham

#Distance of the nucley from grid center
a = 1

#Nuclear charges on centers AB
Za  = 6
Zb = 0

#Set polaization. 1 Unpolarized, 2 Polarized
pol = 2

Nmo = [[2,2], [1,0]]
N   = [[2,2], [2,0]]

optKS = {
        "interaction_type" : "dft",
        "SYM" : False,
        "FRACTIONAL" : True,
        }

#Grid Options
NP = 7 #Number of points per block
NM =  [10,10] #Number of blocks [angular, radial]
L = np.arccosh(15./a) #Maximum radial coordinate value
loc = np.array(range(-4,5)) #Non inclusive on upper bound

#Create and initialize grid object
grid = Psgrid(NP, NM, a, L, loc)
grid.initialize()

#Kohn Sham object
KS = Kohnsham(grid, Za, Zb, pol, Nmo, N, optKS)
KS.scf(optKS)

print(f" Total Energy: {KS.E.E}")