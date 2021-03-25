Prolate Spheroidal Coordinates
==============================

The Prolate Spheroidal Coordinate (SP) system is a two-center orthogonnal coordinate system. It consists of confocal ellipsoids and hyperboloids of revolution. 

Start from the usual Cartesian coordinates and label the foci *a* in the coordinates (0,0,-a) and (0,0,a). The transformation equations are:

.. math:: x = a sinh(u) sin(v) cos(\phi) 
.. math:: y = a sinh(u) sin(v) sin(\phi) 
.. math:: z = a cosh(u) cos(v)

Where u and v label the ellipsoidal and hyperboloidal surfaces, respectively and phi is the azimuth angle. The normal domains of the u and v are:

.. math:: u \in [0, \infty)
.. math:: v \in [0, \pi]

In the case of an axially symmetric system, such is the case of a diatomic molecule, the phi dependance can be treated analytically, and we tehrefore have a two-dimensional problem. 

The reasoning behing choosing this particular coordinate system is that by defining r_1 and r_2 the distances from any given point (x, y, z) to the previously mentioned foci, respectively, we have:

.. math:: r_1 = a (cosh(u) + cos(v))
.. math:: r_2 = a (cosh(u) - cos(v))

And exponentials of the previous functions are transformed into Gaussian functions of u and v in the regions near the nuclei. In other words, there are no exponential nuclear cusps in the (u, v) coordinate system. 

Secondly, although we have previously defined the domains of u and $v$, the functions are well behaved for any real values of them. This allows us to impose usefull constraints on molecular functions at the boundary regions. For example, functions that are of sigma symmetry are even functions of u and v, and same for the odd case (pi). 

In the actual code, it is useful to introduce a new set of coordinates in order to optimize the distribution of the u and the v coordinate surfaces in real space. The transformations are:

.. math:: u = C_1 tanh^{-1}(x_1)
.. math:: v = x_2 + C_2 sin (2x_2)