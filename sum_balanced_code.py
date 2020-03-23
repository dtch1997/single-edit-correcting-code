# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 01:31:03 2020

@author: Daniel Tan
"""


import numpy as np
from qary_string import QaryString
from scipy.special import comb

def binom(N, r):
    return comb(N, r, exact=True)

class CombinatorialBitstringEncoder:
    def __init__(self):
        pass
        
    @staticmethod
    def encode(s, verbose=False):
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
            term = binom(ct, t+1)
            if verbose: print(term, ct, t+1)
            sum += term
        index = sum
        return k, w, index
    
    @staticmethod
    def decode(k, w, index, verbose=False):
        indices = []
        
        # get the maximum index
        for t in range(w, 0, -1):
            ct = 0
            term = binom(ct, t)
            while term <= index and ct < k:
                ct += 1
                term = binom(ct, t)
            ct = ct-1
            term = binom(ct, t)
            if verbose: print(term, ct, t)
            indices.append(ct)
            index = index - binom(ct, t)
        indices = np.array(indices, dtype=int)
        s = np.zeros(k)
        s[indices] = 1
        return s.astype(int)
    
    @staticmethod
    def num_bitstring(k, w):
        # The number of k-length bitstrings with sum w
        return int(binom(k, w))

class SumBalancedCode:
    def __init__(self, k):
        self.k = k
        self._compute_buckets()
        
    def _compute_buckets(self, verbose=False):
        """
        Precompute all sum-pairs (a,b) for which the bitstrings indicate a non-k-
        sum-balanced 4-ary string. 
        
        Note: This code has not been designed to work for q=/=4
        """
        self.sumpair2bucket = {}
        self.bucket2sumpair = []
        self.bucket2size = []

        idx = 0
        k = self.k
        for a in range(k):
            sz_a = CombinatorialBitstringEncoder.num_bitstring(k,a)
            for b in range(k):
                sz_b = CombinatorialBitstringEncoder.num_bitstring(k,b)
                if 2*a + b <= k or 2*a + b >= 2*k:
                    if verbose: print(a, b, idx)
                    self.sumpair2bucket[(a,b)] = idx
                    self.bucket2sumpair.append((a,b))
                    self.bucket2size.append(sz_a * sz_b)
                    idx += 1
                    
        self.bucket2startidx = [0]*len(self.bucket2size)  
        start_idx = 0
        for bucket_idx, bucket_size in enumerate(self.bucket2size):
            self.bucket2startidx[bucket_idx] = start_idx
            start_idx += bucket_size
        return self
    
    @property
    def _num_fwords(self):
        return self.bucket2startidx[-1] + self.bucket2size[-1]
    
    def _fword_to_index(self, word):
        # word: A non-k-sum-balanced word of length k. 
        # qary string (k,) -> binary matrix (k, q) -> bucket (q,) -> bucket index: int
        k, index_in_bucket = CombinatorialBitstringEncoder.encode(word)
        assert k == self.k
        q, bm = word.as_binary_matrix
        sumpair = np.sum(bm, axis=0)
        bucket = self.sumpair2bucket[tuple(sumpair)]
        return self.bucket2startidx[bucket] + index_in_bucket
    
    def _index_to_fword(self, index):
        # index: An integer representing a non-k-sum-balanced word of length k.  
        
        # First locate the bucket
        bucket = 0
        breakpoint()
        while self.bucket2startidx[bucket] < index:
            bucket += 1
        bucket -= 1 # Maximum bucket not exceeding. 
        a,b = self.bucket2sumpair[bucket]
        w = 2*a + b
        index_in_bucket = index - self.bucket2startidx[bucket]
        word = CombinatorialBitstringEncoder.decode(self.k, w, index_in_bucket)
        return word
        
    def encode(self, s):
        """
        Parameter:
            s: An arbitrary QaryString. 
        Return
            x: a k-sum-balanced QaryString
        """
        # Step 1: Append 0
        x = s.concatenate(QaryString(s.q, [0]))
        
        # Step 2: Sequence replacement of all forbidden words
        i = 0 # Index into x
        k = self.k
        
        while i < x.length - k:
            word = x[i:i+k]
            if not word.is_sum_balanced:
                index = self._fword_to_index(word)
                x = x[:i].concatenate([
                    x[i+k:], 
                    QaryString(x.q).fromint(index).pad_to(x.bitlen(self._num_fwords)),
                    QaryString(x.q).fromint(i).pad_to(x.bitlen(s.length)),
                    QaryString(x.q, val=[3])
                ])
            else: 
                i += 1
                
        return x, s.length
    
    def decode(self, x, s_len):
        """
        Parameter:
            x: a k-sum-balanced Qarystring 
        Return
            x: The decoded string
        """
        sentinel = QaryString(x.q, [0])
        # breakpoint()        
        lengths = [x.bitlen(self._num_fwords), x.bitlen(s_len), 1]
        block_len = np.sum(lengths)
        while x[-1] != sentinel:
            index, i, _ = x[-block_len:].split(lengths)
            word = self._index_to_fword(index.asint())
            x = x[:i].concatenate([
                word, 
                x[i:-block_len]
            ])
        x = x[:-1]
        return x
        

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
    code = SumBalancedCode(16)
    for i in range(1000):
        s = QaryString(q=4, val=np.random.randint(low=0, high=4, size=100))
        x, l = code.encode(s)
        s_pred = code.decode(x, l)
        assert s == s_pred