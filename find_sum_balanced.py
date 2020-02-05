# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 00:04:07 2020

@author: Daniel Tan
"""


import numpy as np

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