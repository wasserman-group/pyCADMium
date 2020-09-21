import numpy as np

from CADMium import Psgrid
from CADMium.kohnsham import scf
from CADMium import Kohnsham

#Distance of the nucley from grid center
a = 2.354/2

#Nuclear charges on centers AB
Za  = 6
Zb = 6

#Set polaization. 1 Unpolarized, 2 Polarized
pol = 1

Nmo = [[4], [1]]
N   = [[8], [4]]

optKS = {
        "interaction_type" : "dft",
        "SYM" : True,
        "MaxIter" : 100,
        "ITERLINSOLVE" : False
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

print(f"Total Energy: {KS.E.E}")