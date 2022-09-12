

<p align="center">
<img src="https://raw.githubusercontent.com/wasserman-group/pyCADMium/main/docs/pycadmium_logo_2.png" alt="logo" height=300>
</p>


[![CI](https://img.shields.io/github/workflow/status/wasserman-group/CADMium/CI)](https://github.com/wasserman-group/CADMium/actions?query=workflow%3ACI)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04459/status.svg)](https://doi.org/10.21105/joss.04459)

- Prolate spheroidal coordinates grid-based code that solves a cylindrical problem semi-analytically.  
- Calculations of atoms and diatomic molecules in real space.
- Kohn-Sham DFT calculations using density functional approximations by interfacing with libxc. 
- Density-to-potential inversion calculations. 
- Fragment-based calculations using partition-DFT. 


### Getting started:  
- Clone and install from this repository.
- If installing in Windows, we recommend the use of [WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10).

```
git clone https://github.com/wasserman-group/CADMium.git
cd CADMium
pip install . 
```
pylibxc must be installed as well. 
```
pip install pylibxc2
```
- Alternative, one can install libxc through conda and install pylibxc manually. 
```
conda install -c conda-forge libxc
wget http://www.tddft.org/programs/libxc/down.php?file=5.0.0/libxc-5.0.0.tar.gz # Replace with most recent libxc download link
tar -xf libxc-5.0.0.tar.gz                                                      # Replace with name of downloaded file
cd libxc-5.0.0                                                                  # Replace with name of downloaded file
python setup.py install
```

- Proceed to test your pyCADMium installation. 
Within the pyCADMium/tests/ folder run:
```
pytest .
```
This will run a set of tests dedicated to ensure that your code is fully functional. 
Be wary of warnings, but don't fret on them if present. 

### Community Guideliness
- If any unexpected error occurs, please contact us at: victorandscience@gmail.com.  
- Any comment, suggestion or change is welcomed as an issue or a pull request. 

### Tutorials:
Learn how to use CADMium with these [examples](https://github.com/wasserman-group/CADMium_examples).  
  
### Copyright
Copyright (c) 2020, Wasserman Group  

#### Acknowledgements
*Victor H. Chavez* was supported by a fellowship from The Molecular Sciences Software Institute under NSF grant OAC-1547580.  
Project based on the [MolSSI Cookiecutter](https://github.com/molssi/cookiecutter-cms).  
