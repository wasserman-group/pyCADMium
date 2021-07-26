CADMium
==============================
 *Chemical Atoms in Diatomic Molecules*


[//]: # (Badges)

[![CI](https://img.shields.io/github/workflow/status/wasserman-group/CADMium/CI)](https://github.com/wasserman-group/CADMium/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/wasserman-group/CADMium/branch/master/graph/badge.svg)](https://codecov.io/gh/wasserman-group/CADMium)
<!-- [![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/wasserman-group/CADMium.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/wasserman-group/CADMium/context:python) -->


- Prolate spheroidal coordinates grid-based code that solves a cylindrical problem analitically.  
- Calculations in atoms and diatomic molecules free from basis set incompleteness error. 
- Kohn-Sham DFT calculations using density functional approximations from libxc. 
- Density-to-potential inversion calculations. 
- Fragment-based calculations using Partition-DFT. 


### Getting started:  
- Clone and install from this repository.
- If installing in Windows, we recommend the use of [WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10).

```
git clone https://github.com/wasserman-group/CADMium.git
cd CADMium
pip install . 
```
- Libxc and pylibxc must be installed as well. 
```
conda install -c conda-forge libxc
```
- To communicate libxc with your python site-packages folder:
```
wget http://www.tddft.org/programs/libxc/down.php?file=5.0.0/libxc-5.0.0.tar.gz
cd libxc-5.0.0
tar -xf libxc-5.0.0.tar.gz
python setup.py install
```

- If any unexpected error occurs, please contact us at: gonza445@purdue.edu  

### Tutorials:
Learn how to use CADMium with these [examples](https://github.com/wasserman-group/CADMium_examples).  
  
### Copyright
Copyright (c) 2020, Wasserman Group  

#### Acknowledgements
*Victor H. Chavez* was supported by a fellowship from The Molecular Sciences Software Institute under NSF grant OAC-1547580.  
Project based on the [MolSSI Cookiecutter](https://github.com/molssi/cookiecutter-cms).  
