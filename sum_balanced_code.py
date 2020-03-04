# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 01:31:03 2020

@author: Daniel Tan
"""


import numpy as np
from qary_string import QaryString
from scipy.special import binom

class CombinatorialBitstringEncoder:
    def __init__(self):
        pass
        
    @staticmethod
    def encode(s):
        """
        parameters:
        s: A k-length bitstring. numpy array. 
        
        Return:
            k: length of the string. 
            w, the sum of elements in s, i.e. number of 1's. 
            index, the index of s in the set of k-length bitstrings with sum w. 
        """
        # convert a bitstring to its indices
        # convert the indices to binomial coefficients
        indices = np.nonzero(s)[0]
        w = np.sum(s)
        k = len(s)
        
        sum = 0
        for t, ct in enumerate(indices):
            sum += binom(ct, t+1) 
        index = int(sum)
        return k, w, index
    
    @staticmethod
    def decode(k, w, index):
        indices = []
        
        # get the maximum index
        for t in range(w, 0, -1):
            ct = 0
            while binom(ct, t) <= index:
                ct += 1
            ct = ct-1
            indices.append(ct)
            index = index - binom(ct, t)
        indices = np.array(indices)
        s = np.zeros(k)
        s[indices] = 1
        return s.astype(int)
        

class SumBalancedCode:
    def __init__(self, k):
        self.k = k
        self._compute_buckets()
        
    def _compute_buckets(self):
        """
        Precompute all sum-pairs (a,b) for which the bitstrings indicate a non-k-
        sum-balanced 4-ary string. 
        """
        self.sumpair2bucket = {}
        self.bucket2sumpair = []
        self.bucket2size = {}

        idx = 0
        k = self.k
        for a in range(k):
            sz_a = self._num_bitstring(a)
            for b in range(k):
                sz_b = self.num_bitstring(b)
                if 2*a + b <= k or 2*a + b >= 2*k:
                    self.sumpair2bucket[(a,b)] = idx
                    self.bucket2sumpair[idx] = (a,b)
                    self.bucket2size[idx] = sz_a * sz_b
                    idx += 1
                    
        self.bucket2startidx = []  
        start_idx = 0
        for bucket_idx, bucket_size in self.bucket2size.items():
            self.bucket2startidx[bucket_idx] = start_idx
            start_idx += bucket_size
        return self
    
    def encode(s):
        """
        Parameter:
            s: An arbitrary QaryString. 
        Return
            x: a k-sum-balanced QaryString
        """
        raise Exception("Not implemented yet")
        x = None
        return x
    
    def decode(x):
        """
        Parameter:
            x: a k-sum-balanced Qarystring 
        Return
            s: The decoded string
        """
        raise Exception("Not implemented yet")
        s = None
        return s
        

class SimpleCode:
    """
    Applies a trivial encoding that guarantees k-sum-balancedness
    """
    def __init__(self):
        pass
        
    @staticmethod
    def encode(x):
        """
        x: QaryString
        """
        # Append the inverse of every element. 
        vals = np.zeros([x.length, 2])
        vals[:,0] = x.val
        vals[:,1] = (x.q-1) - x.val
        vals = vals.flatten()
        return QaryString(x.q, np.array(vals))
    
    @staticmethod
    def decode(x):
        # Return elements 0, 2, 4, ... 
        idx = np.arange(start=0, stop=x.length, step=2)
        return x[idx]
    
    
if __name__ == "__main__":
    x = np.array([0,0,1,0,1,1])
    k, w, idx = CombinatorialBitstringEncoder.encode(x)
    s = CombinatorialBitstringEncoder.decode(k,w,idx)
    print(s, x)
    