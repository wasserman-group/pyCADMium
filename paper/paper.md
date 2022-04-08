---
title: 'pyCADMium: Chemical Atoms in Molecules. A prolate spheroidal module for studying atoms and diatomics'  

tags:  
    - partition density functional theory  
    - computational chemistry  
    - density functional theory  
    - embedding methods  

authors:
    name: Victor H. Chávez   
    orcid: 0000-0003-3765-2961    
    affiliations: 1  
    name: Jonathan Nafziger    
    affiliations: 2    
    name: Adam Wasserman    
    orcid: 0000-0002-8037-4453    
    affiliatios: 1, 3  

affiliations:  
    name: Department of Chemistry, Purdue University, West Lafayette, Indiana, USA  
    index: 1  
    name: Gridspace, Los Angeles, California, USA  
    index: 2  
    name: Department of Physics and Astronomy, Purdue Unversity, West Lafayette, Indiana, USA  
    index: 3  

date: 4 April 2022  
bibliography: paper.bib  

---

# Summary 
Diatomic molecules are one of the most exciting systems in quantum chemistry:

- With only a handful of electrons, their computing time is short.
- They highlight the present and ubiquitous problems of modern approximations.
- By defining two fragments with an atom each, diatomics become the canonical system to implement and test embedding methods.

In this work, we introduce ``pyCADMium``, a Python module that uses a prolate spheroidal (PS) coordinate grid to perform computational chemistry calculations on atoms and diatomic molecules accurately. The name is an acronym for "Chemical Atoms in Molecules" that references the idea of chemically significant fragments within a molecule.  
pyCADMium originated in a proprietary programming language but has been rewritten from the ground up as an open-source alternative. The code has been the main driving force behind the development of "Partition Density Functional Theory," [@20CW,@10EW,@17NW,@18JW] a method that uses quantum embedding to lower the cost of a calculation and fixes problems inherent to density functional approximations[@08CY].  

Most practical calculations use a basis-set to represent operators and quantities like potentials, orbitals, and densities [@12SO,@13H]. The lack of an accurate space representation of these and other operators gives rise to the basis-set incompleteness error. This error can be minimized by increasing the number of basis functions used, but it cannot be entirely eliminated in practice. On the other hand, grid-based methods intrinsically allow for an accurate representation [@octopus]. Nevertheless, the number of points required to achieve a significant accuracy can become quickly unmanageable. We wrote ``pyCADMium``, which uses a prolate PS grid to circumvent these issues. The PS grid is a coordinate system formed by revolving an elliptic coordinate plane around the line intersecting the two foci. These planes are formed by ellipses and hyperbolae that share the same focus [@arfken]. Atoms and diatomics are ideally suited for this coordinate system given that each foci can be used to allocate an atom. Additionally, the PS grid is denser around the foci in its cartesian representation, where we expect functions of molecular systems to change rapidly [@17RS].  

Consider a PS grid with foci located at $(-a/2/0)$ and $(a/2,0)$, where $a$ represents the separation distance of a diatomic of interest. Place one (A, ∅) or two atoms (A, B) on each of these foci and specify its charge, number of electrons, and number of atomic and molecular orbitals. Continue by specifying the symmetry for orbitals to calculate. This implies that the user chooses where to allocate each of the avaliable electrons.  

Once the fragments and/or molecule is defined. ``pyCADMium`` can perform:

1. Density functional theory (DFT) calculation. Choose any density functional approximation—up to the generalized gradient approximation (GGA) rung [@03TS]— from the library of exchange correlation functionals libxc [@libxc], and perform a self-consistent Konh-Sham DFT calculation to find the energies, orbitals and density.  

2. Density-to-Potential inversion. Given any density on the PS grid, perfom a numerical inversion [@22SW] to find the multiplicative external potential that reproduces the input density. This problem is ill-posed when the potential is expressed on a basis-set, but it is well-posed when done on a grid [@18JW2].  

3. Partition-DFT. Given a set of fragment densities, perform a numerical inversion to find the non-additve kinetic potential, that when added to the exact non-additve external-hartree-exchange-correlation potential, becomes the exact and unique embedding potential required for the isolated fragments to fully reproduce the molecular density as well as minimizing the sum of fragment energies.  


PS coordinates have proven to excel at accuracy in calculations using atoms and diatomic molecules [@B82,@L09]. Despite these coordinates being used thoroughly in literature, the options for freely-availiable modules that focuses on embedding applications in pretty much non-existent. We have proven in the past that the addition of different embedding methods such as subsystem-DFT can be easily done [@14JW]. We expect that more developers and research take advantage of our code and allow more embedding methods to be implemented in ``pyCADMium``.  

Additionally, a repository of Jupyter Notebooks has been released that includes different examples of the functionalities availiable in the code as well as being able to reproduce some of the previously published results [@].

# Acknowledgements

Victor H. Chávez thanks Dr. Taylor A. Barnes of The Molecular Sciences Software Insitute for the massive help to set the continious integration for the package as well as multiple conversations about the scope and structure of the code.  
This work was supported by the National Science Foundation under GrantNo. CHE-1900301. Victor H. Chavez was supported by a fellowship from The
Molecular Sciences Software Institute under NSF grant OAC-1547580.  
