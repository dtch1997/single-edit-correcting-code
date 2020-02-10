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
    r = np.arange(x.shape[0]) + 1
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
     
def is_sum_balanced(x, q=4):
    k = x.shape[0]
    return (q // 2 - 1) * k < np.sum(x) < (q // 2)*k

def is_k_sum_balanced(x, k, q=4):
    """
    x: numpy array of shape (n,), in a q-ary alphabet. 
    k: Length of window over which we compute sum-balancedness
    """
    n = x.shape[0]
    for i in range(n+1-k):
        window = x[i:i+k]
        if not is_sum_balanced(window, q):
            return False
    return True