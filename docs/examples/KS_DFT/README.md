# 01_ATOMS

- Perform a basic Kohn-Sham DFT calculation on individual atoms up to Neon.
 
- Highlights some of the primary usages of pyCADMium, especially the molecule specification. 

- Some of the options to define a system are self-explanatory. For example, the value `Z_x` assigns the nuclear charge, while pol states whether or not we are performing a spin unrestricted calculation.
 
- Assigning the electronic structure can get a bit finicky. To do this, one must know the system's orbital configuration. That is the number of orbitals for a respective angular momentum. If only one bracket is in the definition, only `s` orbitals are present; If two brackets are present, `p` orbitals are present, etc. Additionally, if the brackets have one component, alpha_i and beta_i orbital are identical; if two values are present, one corresponds to alpha and the other to beta electrons.
 
- If pol = 1, you must assign 2 electrons per orbital. On the other hand, if pol = 2, you must set 1 electron per orbital. 

- Even if there is only one atom in the system, you still have to specify a charge of `Z_x = 0'.

- This notebook helps perform benchmarks for different functionals using a fine grid.

 
