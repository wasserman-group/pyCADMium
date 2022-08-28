---
title: 'pyCADMium: Chemical Atoms in Diatomic Molecules. A prolate spheroidal Python module for embedding calculations'  
tags:  
    - partition density functional theory  
    - computational chemistry  
    - density functional theory  
    - embedding methods  
authors:
    - name: Victor H. Chávez
      orcid: 0000-0003-3765-2961    
      affiliation: 1  
    - name: Jonathan Nafziger    
      affiliation: 2    
    - name: Adam Wasserman    
      orcid: 0000-0002-8037-4453    
      affiliation: "1, 3"  
affiliations:  
    - name: Department of Chemistry, Purdue University, West Lafayette, Indiana, USA  
      index: 1  
    - name: Gridspace, Los Angeles, California, USA  
      index: 2  
    - name: Department of Physics and Astronomy, Purdue University, West Lafayette, Indiana, USA  
      index: 3  
date: 16 Mayo 2022  
bibliography: paper.bib  

---
# Introduction

Diatomic molecules are among the most useful systems to test new ideas in quantum chemistry:

- They highlight the present and ubiquitous problems of modern approximations [@08CY]. 
- By defining two fragments with an atom each, diatomics are ideal to implement and test quantum embedding methods [@20CW].

Most practical calculations use a finite basis set of fixed functional form to represent operators and quantities like potentials, orbitals, and densities [@12SO;@13H]. Such approaches require analytic expressions for the matrix representation of various operators and often exhibit poor basis set convergence.  The basis set incompleteness error can be minimized by increasing the number of basis functions used, but it cannot be entirely eliminated in practice because of the lack of parameters (for instance, exponents, contraction coefficients, etc.) for large basis sets. On the other hand, grid-based methods intrinsically allow for an accurate representation of operators and permit the use of arbitrarily large basis sets [@octopus]. Nevertheless, the number of points required to achieve a significant accuracy can quickly become unmanageable. ``pyCADMium`` presented in this work, uses a prolate spheroidal (PS) grid to circumvent these issues for atoms and diatomic molecules. 

# Summary 

In this work, we introduce ``pyCADMium``, a Python module that uses a PS coordinate grid to accurately perform computational chemistry calculations on systems with cylindrical symmetry. The name is an acronym for ``Chemical Atoms in Diatomic Molecules``. ``pyCADMium`` originated in a proprietary programming language but has been rewritten from the ground up as an open-source alternative. The code has been the main driving force behind the development of "Partition Density Functional Theory," [@20CW;@10EW;@17NW;@18JW] a method that uses quantum embedding to lower the cost of a calculation and fixes problems inherent to density functional approximations [@08CY].  

The PS grid is a coordinate system formed by revolving an elliptic coordinate plane around the line intersecting the two foci. These planes are formed by ellipses and hyperbolae that share the same focus [@arfken]. Atoms and diatomics are ideally suited for this coordinate system given that each foci can be used to allocate an atom. Additionally, the PS grid is denser around the foci in its cartesian representation, where we expect wave functions of molecular systems to change rapidly [@17RS].  

# Functionality

Consider a PS grid with foci located at $(-a/2,0)$ and $(a/2,0)$, where $a$ represents the inter-nuclear separation in a diatomic of interest. Place one (A, ∅) or two atoms (A, B) at each of these foci and specify its charge, number of electrons, and number of atomic and/or molecular orbitals. Continue by specifying the symmetry of orbitals to calculate.

Once the set of fragments and/or molecule is defined ``pyCADMium`` can perform:

1. Density functional theory (DFT) calculation. Choose a density functional approximation up to the generalized gradient approximation (GGA) rung [@03TS] from the library of exchange correlation functionals, libxc [@libxc], and perform a self-consistent Kohn-Sham DFT calculation to find the energies, orbitals and density.  

2. Density-to-Potential inversion. Given any density on the PS grid, perform a numerical inversion [@22SW] to find the multiplicative external potential that reproduces the input density[@03WY]. This problem is ill-posed when the potential is expressed on a basis-set, but it is well-posed when done on a grid [@18JW2].  

3. Partition-DFT. Given a molecule of interest, fragment its external potential and find the set of densities associated with them. The sum of these densities, for most cases, will not be equal to the density of the full system, but we want them to be. Perform a numerical inversion to find the non-additive kinetic potential, that when added to the exact non-additive external-hartree-exchange-correlation potential becomes the unique embedding potential required for the sum of isolated fragments to reproduce the molecular density as well as minimize the sum of fragment energies [@10EW].  

# Statement of Need

PS coordinates have proven to be accurate in calculations using atoms and diatomic molecules [@B82]. Despite these coordinates being used thoroughly in literature, the options for freely-available modules that focus on embedding applications are almost non-existent. One example worth mentioning is ``DARSEC`` [@darsec], a FORTRAN-based code that uses the same prolate spheroidal methodology in the context of the optimized effective potential equation. Their use of the coordinates and basic operators, such as derivatives and integration, is analogous to our approach. A key difference is our use of `libxc` to obtain the functional derivatives needed to obtain exchange-correlation potentials. Moreover, their codebase is not openly available. 

Our code was designed from the ground up with embedding calculations in mind. As described in [@14NW], it can be used to perform and quickly develop other embedding calculations. Last but not least, a repository of Jupyter Notebooks is available on [Github](https://github.com/wasserman-group/CADMium_examples) that includes various examples of the functionalities available in the code. 

# Acknowledgements

Victor H. Chávez thanks Dr. Taylor A. Barnes of The Molecular Sciences Software Institute for the massive help to set the continuous integration for the package as well as multiple conversations about the scope and structure of the code.  
This work was supported by the National Science Foundation under Grant No. CHE-1900301. Victor H. Chávez was supported by a fellowship from The Molecular Sciences Software Institute under NSF grant OAC-1547580.  

# References

