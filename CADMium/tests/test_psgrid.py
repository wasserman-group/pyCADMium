"""

test_psgrid.py
Tests attributes and methods for an instance of psgrid

"""

import pytest
from CADMium.psgrid import Psgrid
import numpy as np

@pytest.fixture()
def grid_test():
    
    #Distance of the nuclei from grid center
    a = 0.5
    #Grid options
    NP = 7
    #Number of blocks
    NM = [4,4]
    #Maximum radial coordinate
    L = np.arcsinh(13.0/a)
    loc = np.array(range(-4,5))
    #Initialize
    grid = Psgrid(NP, NM, a, L, loc)
    grid.initialize()

    return grid


def test_odd_angular_derivatives_1D(grid_test):
    """
    Angular Derivatives
    xa (angular coordinate) goes from 0 to pi. 
    Sin has odd boundary conditions at each end point so we use oDa1
    """

    numerical = grid_test.oDa1 @ np.sin(grid_test.xa)
    analytic = np.cos(grid_test.xa) #Derivative of sin is cos
    assert(np.isclose(np.abs(numerical).all(), np.abs(analytic).all()))

def test_odd_angular_derivatives_2D(grid_test):
    """
    Angular Derivatives
    xa (angular coordinate) goes from 0 to pi. 
    Sin has odd boundary conditions at each end point so we use oDa2
    """

    numerical = grid_test.oDa2 @ np.sin(grid_test.xa)
    analytic = 1-.0 * np.sin(grid_test.xa) #Derivative of sin is cos
    assert(np.isclose(np.abs(numerical).all(), np.abs(analytic).all()))

def test_even_angular_derivatives_1D(grid_test):
    """
    Angular Derivatives
    xa (angular coordinate) goes from 0 to pi. 
    Sin has odd boundary conditions at each end point so we use eDa1
    """

    numerical = grid_test.eDa1 @ np.cos(grid_test.xa)
    analytic = -1.0 * np.sin(grid_test.xa) #Derivative of sin is cos
    assert(np.isclose(np.abs(numerical).all(), np.abs(analytic).all()))

def test_even_angular_derivatives_2D(grid_test):
    """
    Angular Derivatives
    xa (angular coordinate) goes from 0 to pi. 
    Sin has odd boundary conditions at each end point so we use eDa2
    """

    numerical = grid_test.oDa2 @ np.cos(grid_test.xa)
    analytic = -1.0 * np.cos(grid_test.xa) #Derivative of sin is cos
    assert(np.isclose(np.abs(numerical).all(), np.abs(analytic).all()))
    