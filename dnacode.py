# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 23:28:57 2020

@author: Daniel Tan
"""

import numpy as np

def syndrome(x):
    assert len(x.shape) == 1
    r = np.arange(x.shape[0]) + 1
    return np.sum(x * r)

def signature(x):
    assert len(x.shape) == 1
    return x[1:] >= x[:-1]

class DNACode:
    def __init__(self):
        pass
    
    def encode(self, x):
        """
        Parameters
        ----------
        x : numpy array of shape (n-1,)
            Integers in a q-ary alphabet, ranging from 0 to q-1. 

        Returns
        -------
        x_enc : numpy array of shape (n,)
                Integers in a q-ary alphabet. 

        """
        x_enc = None
        ###
        y = None # Get the the corresponding codeword in Bal_k(x). 
        
        ###
        return x_enc
    
    def decode(self, x_enc):
        x = None
        return x