"""
CADMium
A short description of the project.
"""

# Add imports here
from .psgrid import Psgrid
from .inverter import Inverter
from .pssolver import Pssolver
from .partition import Partition
from .kohnsham import Kohnsham

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
