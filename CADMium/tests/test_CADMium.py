"""
Unit and regression test for the CADMium package.
"""

# Import package, test suite, and other packages as needed
import CADMium
import pytest
import sys

def test_CADMium_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "CADMium" in sys.modules
