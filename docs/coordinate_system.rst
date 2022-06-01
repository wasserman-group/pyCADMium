The Coordinate System
=====================

The PS grid is a two-center orthogonal coordinate system formed by revolving an elliptic coordinate plane around the line intersecting two foci. These planes are formed by ellipses and hyperbolae that share the same focus \cite{arfken}. Assume that we place the foci in the Cartesian coordinates  :math:`(0,0,-a)` and :math:`(0,0,a)`. We can transform from one coordinate system to the other with

.. math:: x = a \cdot sinh(\mu) \cdot sin(\nu) \cdot cos(\phi) 
.. math:: x = a \cdot sinh(\mu) \cdot sin(\nu) \cdot sin(\phi) 
.. math:: x = a \cdot cosh(\mu) \cdot cos(\nu)

Where :math:`a` is the distance between the foci, and :math:`\mu` and :math:`\nu` correspond to the ellipsoidal and hyperboloid surfaces. In the case of an axially symmetric system, the molecular orbitals of interest are independent of the azimuth angle according to Equations (\ref{coordinate_transformations}). The :math:`\phi` coordinate can be treated analytically and our problem now becomes two-dimensional. 

To allow a more accurate description of the systems around the nuclei, the coordinates can be mapped into a new set of 

.. \newcommand{\sate}{\mathlarger{\mathfrak{e}}}
.. \newcommand{\sath}{\mathlarger{\mathfrak{h}}}

.. math:: \mathfrak{e} = cosh^{-1}(\mu), \quad \mathfrak{e} \in [0, \infty)
.. math:: \mathfrak{h} = cos^{-1}(\nu), \quad  \mathfrak{h} \in [0,\pi]

In the new coordinates, the 'radial' part of the laplacian is 

.. math::  \nabla^2 =  \frac{4}{a^2 (\mu^2 - \nu^2)} \Bigg\{ \frac{\partial{^2}}{\partial{\mathfrak{e}^2}}
                                           \frac{\mu}{\sqrt{\mu^2 - 1}} \frac{\partial}{\partial \mathfrak{e}}
                                           +\frac{\partial^2}{\partial \mathfrak{h}^2}
                                           +\frac{\nu}{\sqrt{1-\nu^2}} \frac{\partial}{\partial \mathfrak{h}}
                                           - \mathbf{m_a^2} \cdot \Bigg(  \frac{1}{\mu - 1} + \frac{1}{1 - \nu^2}  \Bigg)
                                           \Bigg\}


where :math:`\mathbf{m_a^2}` is an integer and it specifies the rotation symmetry of the orbitals \cite{96KS}:  \\

.. math:: m_a = 0      \rightarrow \sigma \quad orbitals
.. math:: m_a = \pm 1  \rightarrow \pi    \quad orbitals
.. math:: m_a = \pm 2  \rightarrow \delta \quad orbitals
.. math:: m_a = \pm 3  \rightarrow \phi   \quad orbitals


Orbitals that have a higher symmetry than :math:`\phi` are not relevant for diatomic molecules at the DFT level of theory. Orbitals that have the same radial part and the same :math:`m=\pm m_a` belong to the same shell. 


Symmetry Considerations
-----------------------

These coordinates have a convenient correspondence with the distances between any point at either of the foci of the coordinate system.

.. math:: r_1 = 2a ( cosh(\mathfrak{e}) + cos(\mathfrak{h}) ) = a (\mu + \nu)
.. math:: r_2 = 2a ( cosh(\mathfrak{e}) - cos(\mathfrak{h}) ) = a (\mu - \nu)  

where :math:`r_1` and :math:`r_2` represent the distances of any point :math:`(x,y,z)` to the foci. Additionally, the scales factors and the volume element for a PS grid are

.. math::    h_1 = h_2 = (r_1 \cdot r_2)^{1/2}
.. math::    h_3 = a \cdot sinh(\epsilon) \cdot sin(\hbar)

Further discussions on the scale factors can be found on the curvilinear coordinate system section of Reference \cite{arfken}.
Equations \ref{eq:rinPS} reveal a lot of interesting properties of the transformation in Equations \ref{eq:sat} and \ref{eq:sath}. If the sign of :math:`\mu` or :math:`nu` is changed, then the point :math:`(x,y,z) \rightarrow (-x,-y,z)`. Thus, a rotation by :math:`\pi` radians leaves orbitals with even :math:`m` unchanged  but reverses the sign of the orbitals with odd m-values:

    - Functions with :math:`\sigma, \delta, ...` symmetry are even functions of :math:`(\mu, \nu)`
    - Functions with :math:`\pi, \phi, ...` symmetry are odd functions of :math:`(\mu, \nu)`

Mathematically, this is expressed as

.. math::    f(\mu, \nu)       = (-1)^m f(\mu, -\nu)
.. math::    f(\mu, \nu)       = (-1)^m f(-\mu, \nu)
.. math::    f(\pi + \mu, \nu) = (-1)^m f(\pi-\mu, \nu)

Due to these properties, we are able to use the central difference formula in any function near the boundary line. These relations are also used to adjust :math:`f` along :math:`(0,\mu), (\pi, \mu)` and :math:`(\mu, 0)` boundary lines for :math:`\sigma` functions. Any function that has a higher :math:`\sigma` symmetry vanish at these boundary lines.  


Discretization
--------------


To discretize our system, we build a two-dimensional mesh with two evenly spaced one-dimensional grids. Since this system is orthogonal, the mesh is a product of two independent one-dimensional meshes that generate a rectangle of :math:`N_a \times N_r` points, representing the number of angular and radial points. With appropriate boundary conditions, only half of the rectangle is required for homonuclear diatomics, otherwise the whole rectangle is constructed. After having defined our space, we require our differentiation and integration operations. Integration is performed through

.. math:: \int f(\mathbf{r}) \quad d\mathbf{r} \simeq \sum_{ij}^{N_1} w_{ij} \quad f_{ij}


To calculate the weights proceed as stated in the Appendix \ref{ap:Newton}. The first and second derivatives of the Laplacian defined on Equation \ref{eq:laplacian} are approximated by finite different expression from the central difference formula. 

.. math::     \frac{du}{dx}_j = \frac{1}{dx} \cdot (c1 \cdot u_{[i-j]} + c2 \cdot u_{[j]} + c3 \cdot u_{[j+1]})

where :math:`\{c_i\}` are the finite difference coefficients.  An additional requirements for these operators is that they must have the appropriately boundary conditions for each symmetry class, therefore a distinct matrix is needed for each of the different orbitals symmetries. 

Although we have computed the integration and differentiation operators, we have not discuss how to produce the potentials that are usually used in a KS-DFT calculation. Each of them needs to be addressed individually. 

    - The **external** potential is found by substitution of the transformation equations.
    .. math:: v_{ext}(\mathbf{r}) = - \frac{ Z (cosh(\mathfrak{e}) - cos(\mathfrak{h}))}{ a (cosh^2(\mathfrak{e}) - cos^2(\mathfrak{h})) }

    - The **Hartree** or the Coulomb potential arises from the total electron density. Although we know its explicit form as an integral involving the electronic density, it is more commonly known to be determined by solving the Poisson equation 
    .. math:: \nabla^2 v_{H}(\mathbf{r}) = -4\pi n(\mathbf{r}) 
    The Poisson equation is not solved by direct substitution of the matrices that we produced in the previous section. Instead they are solved using an LU decomposed Laplacian followed by setting up the boundary conditions as suggested in Reference \cite{96KS}.

    - The **exchange-correlation** potential is found through the library of exchange-correlation functionals *libxc* \cite{libxc}. The library does not care about the coordinate system of the density, it only requires as input the value of the density and the gradient at the total discretized points. Currently no higher rungs than GGA are available in the current version of \cad.