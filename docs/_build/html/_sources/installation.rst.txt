Installation
============

- Clone and install from `this <https://github.com/wasserman-group/pyCADMium>`_ repository.  
- If installing in Windows, we recommend the use of `WSL <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_.  

.. code-block:: bash

    git clone https://github.com/wasserman-group/CADMium.git
    cd CADMium
    pip install . 

pylibxc must be installed as well.  

.. code-block:: bash

    pip install pylibxc2

- Alternative, one can install libxc through conda and install pylibxc manually. 

.. code-block:: bash

    conda install -c conda-forge libxc
    wget http://www.tddft.org/programs/libxc/down.php?file=5.0.0/libxc-5.0.0.tar.gz # Replace with most recent libxc download link
    tar -xf libxc-5.0.0.tar.gz                                                      # Replace with name of downloaded file
    cd libxc-5.0.0                                                                  # Replace with name of downloaded file
    python setup.py install

