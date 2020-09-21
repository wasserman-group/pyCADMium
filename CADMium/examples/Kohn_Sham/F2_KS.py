import numpy as np

from CADMium import Psgrid
from CADMium.kohnsham import scf
from CADMium import Kohnsham

#Distance of the nucley from grid center
a = 2.615/2

#Nuclear charges on centers AB
Za  = 9
Zb = 9

#Set polaization. 1 Unpolarized, 2 Polarized
pol = 1

Nmo = [[5], [2]]
N   = [[10], [8]]

optKS = {
        "interaction_type" : "dft",
        "SYM" : True,
        "xc_family" : 'lda', 
        "x_func_id" : 1,
        "c_func_id" : 12,    
        }
#Grid Options
NP = 7 #Number of points per block
NM =  [6,6] #Number of blocks [angular, radial]
L = np.arccosh(15./a) #Maximum radial coordinate value
loc = np.array(range(-4,5)) #Non inclusive on upper bound

#Create and initialize grid object
grid = Psgrid(NP, NM, a, L, loc)
grid.initialize()

#Kohn Sham object
KS = Kohnsham(grid, Za, Zb, pol, Nmo, N, optKS)
KS.scf(optKS)

print(f"Total Energy: {KS.E.E}")