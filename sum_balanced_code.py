# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 01:31:03 2020

@author: Daniel Tan
"""


import numpy as np
from qary_string import QaryString

class SumBalancedCode:
    """
    Applies a trivial encoding that guarantees k-sum-balancedness
    """
    def __init__(self, q: int):
        self.q = q
        
    def encode(self, x):
        """
        x: QaryString
        """
        # Append the inverse of every element. 
        vals = np.zeros([x.length, 2])
        vals[:,0] = x.val
        vals[:,1] = (self.q-1) - x.val
        vals = vals.flatten()
        return QaryString(self.q, np.array(vals))
    
    def decode(self, x):
        # Return elements 0, 2, 4, ... 
        idx = np.arange(start=0, stop=x.length, step=2)
        return x[idx]
            
if __name__ == "__main__":
    x = QaryString(q=4, val=[0,1,2,3])
    code = SumBalancedCode(x.q)
    x_enc = code.encode(x)
    print(x_enc)
    x_pred = code.decode(x_enc)
    assert x == x_pred
    
    