# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 00:10:55 2020

@author: Daniel Tan
"""

import numpy as np

def syndrome(x):
    """
    Computes sum( i*x_i) 
    x: numpy array 
    n: Modulus of the syndrom. 
    """
    assert len(x.shape) == 1
    r = np.arange(x.shape[0])
    return np.sum(x * r)

def signature(x):
    assert len(x.shape) == 1
    return x[1:] >= x[:-1]

def parity_check(x):
    return np.sum(x) % 2

def multiplicative_inverse(a, n):
    """
    Return b such that ab = 1 (mod n)
    """
    for b in range(n):
        if a * b % n == 1: return b
    # Should never get here
    raise Exception()
     