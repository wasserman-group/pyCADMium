Getting Started
===============

In this next section, we briefly go over each of the steps that one needs to do to use the many features in *pyCADMium*. We will start with defining the PS grid. To define a grid within *pyCADMium*, we need to specify the following parameters

.. code-block:: python

    from pycadmium import Psgrid

    a  = 2.615/2                    # Distance between the two foci / 2 (a.u.)  
    NP = 7                          # Number of points per block  
    NM = [6,6]                      # Number of blocks [angular, radial]  
    L = np.arccosh(max_rad/a)       # Maximal radial coordinate value  
    loc = np.array(range(-4 , 5))   # Stencil for finite difference  

    grid = Psgrid(NP, NM, a, L, loc)   

where :math:`a` is half the distance between the two foci. :math:`NP` refers to the number of points in a block, :math:`NM` is the number of points per block, given separately for the 'angular' and 'radial' one dimensional grids to define the PS rectangle. 
The blocking of the points does not affect the calculation numerically, it only multiplies each of the values of :math:`NM` so that the total number of points in the grid equals :math:`(NP-1) * (NM[0]) * (NP-1) * (NM[1])`. 
We defined the size of the box as the maximal radial coordinate value. This quantity is modified according to the *coordinate system* section . The value of :math:`L` must be contained within :math:`(-1,1)`. `loc` is required to select the coefficients for the finite difference approximation. 
Lastly, we provide all of the defined elements to generate a ``PSgrid`` object and we proceed to initialize it.   


---

Molecules in the PSgrid
-----------------------


Let us focus on the first row homonuclear diatomics B :math:`_2`, C :math:`_2`, N :math:`_2`, O :math:`_2` and F :math:`_2`. To use *pyCADMium* we stated that one needs to know about the orbital configuration of the diatomics. The next table shows the orbital configuration for these and other systems. 
 
| **Molecule** | **Orbital Configuration** | **State** |

| B2 | :math:`1\sigma_g^2 \quad 1\sigma_u^2 \quad 2\sigma_g^2 \quad 2\sigma_u^2 \quad 1\pi_u^2` | :math:`^3\Sigma_g^-` |

| C2 | :math:`1\sigma_g^2 \quad 1\sigma_u^2 \quad 2\sigma_g^2 \quad 2\sigma_u^2 \quad 1\pi_u^4` | :math:`^1\Sigma_g^+` |

| N2 | :math:`1\sigma_g^2 \quad 1\sigma_u^2 \quad 2\sigma_g^2 \quad 2\sigma_u^2 \quad 3\sigma_u^2 \quad 1\pi_u^4` | :math:`^1\Sigma_g^+` |

| O2 | :math:`1\sigma_g^2 \quad 1\sigma_u^2 \quad 2\sigma_g^2 \quad 2\sigma_u^2 \quad 3\sigma_g^2 \quad 1\pi_u^4 \quad 1\pi_g^2` | :math:`^3\Sigma_g^-` |

| F2 | :math:`1\sigma_g^2 \quad 1\sigma_u^2 \quad 2\sigma_g^2 \quad 2\sigma_u^2 \quad 3\sigma_g^2 \quad 1\pi_u^4 \quad 1\pi_g^4` | :math:`^1\Sigma_g^+` |  



Let us use the F :math:`_2` molecule as the example of the chapter. To define the geometry in *pyCADMium* we need to provide the following 

.. code-block:: python

    a   = 2.615/2      # Separation distance / 2 (a.u.)
    Za  = 9.0          # Fluorine Atom 1 Nuclear Charne
    Zb  = 9.0          # Fluorine Atom 2 Nuclear Charge
    pol = 1            # Set unpolarized system
    Nmo = [[5],  [2]]  # Number of Molecular Orbitals
    N   = [[10], [8]]  # Number of electrons

In the previous code, the :math:`a` parameter is both the separation distance of the nuclei and the foci of the PS grid. We proceed to define the nuclear charges as well as the polarization of the system. Setting ``pol = 2`` would explicitly compute the :math:`\alpha` and :math:`\beta` components of the density. 
Finally, the molecular configuration along with its symmetries are specied in the next two lines. ``Nmo`` requires the ``total`` number of orbitals per symmetry. 
As we can appreciate from the orbital configuration table, this system has a total of 5 :math:`\sigma` orbitals and 2 :math:`\pi` orbitals, where we allocate 10 and 8 electrons, respectively. 
Nothing prevents us from plugging in different values of the orbitals and the electrons. Thus one must be mindful of the desired configuration. 

---


Kohn-Sham calculation
---------------------

We are ready to finally complete a calculation using *pyCADMium*. The first calculation that we are going to do is to obtain the LDA energy of the F_2 molecule. 
To do so, we need to supply the previous geometry to a ``Kohnsham`` object. Besides the geometry, we can supply an additional dictionary with several options to control how our calculation behaves. 
The simplest example can be constructed as 

.. code-block:: python

    from pycadmium import Kohnsham

    optKS = { # Options for the KS calculation
            'interaction_type' : 'dft',
            'sym'              : True,
            }
    optSCF = { # Options for the SCF cycle
            'maxiter' : 100,
            'e_tol'   : 1e-9,
             }
    
    KS = KohnSham(grid, Za, Zb, pol, Nmo, N, optKS)
    KS.scf(optSCF)


where we requested a KS-DFT calculation. If no functional is provided, by default *pyCADMium* uses the LDA. 
The ``sym`` option specifies that only half of the PS plane needs to be constructed due to the symmetry of the F_2 molecule. 
The next set of options controls the SCF procedure behaviour. Here we specified that we want a maximum of 100 iterations and that the procedure should stop if the differences in energy are less than 1e-9.
We proceed by defining the ``Kohnsham`` object and start the SCF procedure. The results of the calculations such as the energies and potentials is available as properties of the ``Kohnsham`` object.


 P-DFT calculation
 -----------------


One of the highlights of *pyCADMium* is its ability to treat each atom as a fragment. This feature allows it to be used for developing  embedding methods. In this section, we will briefly discuss the algorithm for a P-DFT calculation.
Consider a set of two fragments :math:`\{ n_1, n_2 \}`, each having :math:`\{N_1, N_2\}` electrons respectively. We are interested in recreating the results of a KS-DFT calculation of the full system using only the information of the fragments. 
The number of electrons in each fragment must add up to :math:`N_m`, the number of electrons in the full molecule. Additionally, we seek to minimize the sum of the energies of each fragments. Mathematically, this can be done through the unconstrained minimization of 

.. math::    G[n] = \min_{ \{ n_{\alpha} \} } \Bigg\{ E_f[\{ n_{\alpha} \}] + \int d\mathbf{r} \cdot v_p(\mathbf{r}) \cdot (n_f(\mathbf{r}) - n(\mathbf{r}))  - \sum_{\alpha}\mu_{\alpha} \Bigg( \int d\mathbf{r}  \cdot n_{\alpha}(\mathbf{r}) - N_{\alpha} \Bigg) \Bigg\}

Where each :math:`E_f[n_1]` corresponds to the each of the fragment energies. :math:`E_f[n_1] = T_s[n_1] + E_{\rm Hxc}[n_1] + \int v_{\alpha}(\mathbf{r}) n_1(\mathbf{r}) d\mathbf{r}` and the two other terms are the Lagrange multipliers: the partition potential :math:`v_p(\mathbf{r})` and the chemical potential :math:`\mu_{\alpha}`. 
To minimize, perform a functional derivative with respect to to the density :math:`n_1`

.. math:: \frac{\partial G}{ \partial n_1} = 0
.. math:: \frac{\partial E[\{ n_{\alpha} \}]}{\partial n_1 } + v_p(\mathbf{r}) = \mu_{\alpha}

To solve the above equation, we use a set of KS systems as fragments, each of them in the influence of the multiplicative potential :math:`v_{\alpha}(\mathbf{r}) + v_p(\mathbf{r})`. With the exact partition potential, we can solve the Kohn-Sham equations for the P-DFT problem

.. math:: \{ -\frac{1}{2} \nabla^2 + v_{\alpha,eff}[\{ n_{\alpha} \}](\mathbf{r}) + v_p[\{ n_{\alpha} \}](\mathbf{r}) \} \phi_{i, \alpha}(\mathbf{r}) = \varepsilon_{i, \alpha} \phi_{i, \alpha}(\mathbf{r})

Since the potentials inside the brackets depend on the fragment densities—the quantity we are looking for— we solve these equations self-consistently. At convergence, we obtain a set of fragment orbitals that we use to build the fragment densities as :math:`n_1 = \sum_{i}^{N_f} | \phi_{i, \alpha}(\mathbf{r}) |^2` .

In *pyCADMium* we make use of the class ``Partition`` to define an embedded calculation. This class allows us to define the properties of each fragments as well as the properties of the full molecular system that we ought to solve. To define a ``Partition`` object, follow the F_2 example that we have been studying 

.. code-block:: python

    from pycadmium import Partition
    
    optPartition = { # Options for the object Partition
                     'kinetic_part_type' : 'inversion',
                     'ab_sym'            :  True,
                     'ens_spin_sym'      :  True
                     'isolated'          :  True}}
    
    part = Partition(grid, Za, Zb, pol, MOa, Na, nua, MOb, Nb, nub, optPartition)
    part.scf() 

In the previous example, we first set up the options for our object using a dictionary. By setting ``isolated = True``, we won't be doing a P-DFT calculation just yet. 
Instead, we want to generate an initial guess for the calculation by simply adding the components of the individual fragments. Since the fragments are have no interaction between each other there is no partition potential, and thus no partition energy. 

---

To find the exact partition potential, we make use of numerical inversions. This allows us to obtain the exact Kohn-Sham potential for the sum of fragment densities (This quantity is **not obtained** in a ``forward`` manner, and thus, we don't have access to the Kohn-Sham potential that produces such a density), as well as the kinetic potential of the sum of fragment densities needed to compute the partition potential.

In *pyCADMium*, we make use of the class ``Inverter`` to access a handful of methods to invert a density. For brevity, we will only discuss the *orbital invert* method that has shown to be the most robust and reliable. 
The numerical inversion is essentially a direct search for orbitals that satisfy the KS equations while satisfying some constraints at each point of the two-dimensional space. The conditions form a set of residuals:

    1. Orbitals solve the KS Equations :math:`\rightarrow res_{ij}^{KS} = \{-1/2 \nabla^2 \phi_i\}_j + v_{s,j}\phi_{i,j} - \varepsilon_{i} \phi_{i,j}`
    2. Orbitals are normalized :math:`\rightarrow res_{ij}^{N} = \sum_j \mid \phi_{i,j} \mid^2 - 1`
    3. Sum of the squares of each orbitals equals the target density :math:`\rightarrow res_{ij}^{n} = \sum_j \mid \phi_{i,j} \mid^2 - n_j`

where the :math:`j`` index runs over grid points and the :math:`i`` point runs through the orbitals. Combine all residuals to form a vector function that takes the orbitals, energies and effective potential as its argument. 
If we find the root of this function, we find the effective potential that reproduces the given density. The Jacobian of the residual function happens to be a sparse square matrix that can be treated with the Newton-Raphson minimization to find which vector minimizes the residuals. 
Additionally, the HOMO eigenvalue is fixed to be zero to avoid shifting the potential, and the HOMO normalization constraint is removed allowing it to be satisfied by the overall density constraint. 
These numerical inversion are used at each step in the P-DFT scf procedure to determine the functional derivatives of :math:`ts_{nad}`. 

Remember that we already constructed an initial guess for our P-DFT calculation. Then, in order to complete it alongside with the inversion procedure, we write  

.. code-block:: python

    from pycadmium import Partition
    from pycadmium import Pssovler, Inverter

    optInverter = { # Options for the object Inverter
                    'invert_type'     : 'orbitalinvert',
                   }
    
    # Inverter Requires a PSsolver Object
    mol_solver = Pssolver(grid, MOm, Nm, {}) 
    part.inverter = Inverter(grid, mol_solver, optInverter)

    # We now want the fragments to be under the presence of the partition potential
    part.optPartition.isolated   = False     
    # 'continuing' ensures we use the generated initial guess
    part.scf({"continuing" : True})          


We will be using the ``orbitalinvert``, this can be set in the options for the ``Inverter``. A ``Psolver`` object is needed. This object simply stores the results from a calculation. 
Each ``Kohnsham`` object defines one per orbital. We define the ``Inverter`` *as an attribute of the partition*. We continue by specifying the ``Partition`` object that want to perform a P-DFT calculation without isolated fragments and we trigger the SCF calculation.

In this section, we do not discuss any results. To find a set of examples for different use cases of *pyCADMium*, please visit ``github.com/wasserman-group/pyCADMium_examples``. Where the results for this and many other calculations are found. 
