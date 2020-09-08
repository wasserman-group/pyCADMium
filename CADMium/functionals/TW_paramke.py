"""
TW_paramke.py
"""

def TW_paramke(s, k):

    F  = 1 + k[0] - k[0] / (1 + k[1] / k[0] @ s ** 2 )

    return F