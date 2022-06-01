"""

finite_difference_coefficients.py
TT.py
Determine the coefficients to a finite differencing scheme

"""

import numpy as np
import math

def finite_difference_coefficients(x, deriv):
    """
    Find the coefficients to a finite differencing scheme. 
    i.e. find the coefficients [c1, c2, c3] that are appropriate for the 2nd order
    centered difference approximation to the first derivative:

    (du/dx)[j] = 1/dx * (c1*u[i-j] + c2*u[j] + c3*u[j+1]) + error terms

    Parameters
    ----------
    x: numpy array:
        Vector containing relative distances to sampling points 
    deriv: int
        Order of the derivative to approximate

    Returns
    -------
    coefficients : numpy array
        Coefficients requested
    coefficient_table : numpy array
        Output matrix showing the Taylor table coefficients
    error_order : lowest order of higher order terms 
                  (order of accuracy (e.g. O(dx**2) is 2nd order))
    """

    N = len(x)
    coefficient_table = np.zeros((N, N + 7))

    #Do for each coefficient of the Taylor Expansion
    #Expand out 6 extra terms for higher accuracies at higher derivatives
    for j in range(N):
        for n in range(N + 7):
            coefficient_table[j, n] = -1.0 / math.factorial(n) * x[j]**n

    #Allocate derivative
    d = np.zeros((N, 1))
    d[deriv] = -1

    #Calculate the coefficients
    c = np.linalg.lstsq(coefficient_table[:, 0:N].T, d, rcond=-1)[0]

    #Set any small value to zeros
    for i in range(len(c)):
        if np.abs(c[i]) < 1e-10:
            c[i] = 0
    
    """
    Determine the order of accuracy
    Multiply the Taylor table of coefficients by the finite difference 
    coefficients and add up all like terms. 
    Terms that don't add to zero are the objective derivative and 
    the residual higher order terms
    """
    
    a = (c.T @ coefficient_table)[0]

    #Values closer to zero are zero
    for i in range(len(a)):
        if np.abs(a[i]) < 1e-5:
            a[i] = 0
            
    """
    Location of first nonzero term is the order of derivative you are approximating. 
    Second nonzero term tells the lowest order of the higher order terms
    """
    
    position = np.argwhere(a)[:2]
    # error_order = position[1][0] - deriv
    # coefficient_table = coefficient_table[:, 0:N]  
    
    return c.T[0]