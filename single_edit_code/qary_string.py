# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 00:55:22 2020

@author: Daniel Tan
"""

import numpy as np
from .util import syndrome, signature, parity_check, is_sum_balanced

class QaryString:
    
    def __init__(self, q : int = 4, 
                       val = np.zeros(shape=0, dtype=np.int8)):
        self.q = q
        if type(val) == int or type(val) == np.int8:
            val = [val]
        if type(val) == list:
            val = np.array(val)
        self.val = val.astype(np.int8)
        
    """
    Core methods
    """
        
    @property
    def length(self):
        return self.val.shape[0]
        
    def concatenate(self, qstrs):
        val = np.copy(self.val)
        for qstr in qstrs:
            assert self.q == qstr.q
            val = np.concatenate([val, qstr.val])
        return QaryString(self.q, val)
    
    def asint(self):
        sum = 0
        for i in range(len(self)):
            sum = sum * self.q
            sum += int(self.val[i])
        return sum
    
    def fromint(self, n: int):
        vals = []
        while n > 0:
            vals.append(n % 4)
            n = n // 4
        vals = vals[::-1]
        return QaryString(self.q, val=np.array(vals))
    

    
    """
    Convenience methods
    """
    
    def pad_to(self, length, val=0, at="front"):
        if length < self.length:
            return self
        r = QaryString(self.q, val=np.zeros(length - self.length))
        if at == "front":
            return r.concatenate([self])
        elif at == "back":
            return self.concatenate([r])
        else:
            raise Exception("Invalid value provided for 'at'")
    
    def bitlen(self, n: int):
        # How many symbols in q-ary are required to represent n. 
        # Approximately log_q(n)
        # Need to calculate manually due to overflow issues
        b = 1
        bpow = self.q ** b
        while bpow < n:
            b += 1
            bpow *= self.q
        return b
    
    def split(self, lengths):
        assert np.sum(np.array(lengths)) == self.length
        parts = []
        idx = 0
        for i, l in enumerate(lengths):
            assert l >= 0
            parts.append(self[idx: idx + l])
            idx = idx + l
        return tuple(parts)
    
    def randomize(self, n):
        return QaryString(self.q, np.random.randint(low=0, high=self.q, size=n))
        
    def __len__(self):
        return self.length
    
    def __str__(self):
        return f"q = {self.q}, val = {self.val}"
    
    def __getitem__(self, key):
        return QaryString(self.q, self.val.__getitem__(key))
    
    def __setitem__(self, key, item):
        self.val.__setitem__(key, item)
        
    def __delitem__(self, item):
        self.val.__delitem__(item)
        
    def __eq__(self, other):
        if type(other) == int:
            return self.length == 1 and self.val[0] == other
        return (self.q == other.q and np.all(self.val == other.val))
    
    """
    Methods for performing the encoding
    """
    def mutate(self, mtype="random"):
        pos = np.random.randint(low=0, high=self.length)
        val = np.copy(self.val)
        
        mtypes = ["substitute", "insert", "delete"]
        if mtype == "random":
            mtype = np.random.choice(mtypes)
        if mtype == "substitute":
            allowed_symbols = [i for i in range(self.q)]
            allowed_symbols.remove(val[pos])
            val[pos] = np.random.choice(allowed_symbols)
            symbol = val[pos]
        elif mtype == "insert":
            symbol = np.random.randint(low=0, high=self.q)
            val = np.insert(val, pos, symbol)
        elif mtype == "delete":
            symbol = val[pos]
            val = np.delete(val, pos)            
        return QaryString(self.q, val), mtype, pos, symbol
    
    def _insert(self, pos, symbol, idx_of_pos=0):
        """
        idx_of_pos: If 1, pos is 1-indexed. If 0, pos is 0-indexed.  
        """
        val = np.insert(self.val, pos-idx_of_pos, symbol)
        return QaryString(self.q, val)
    
    def _delete(self, pos, idx_of_pos=0):
        val = np.delete(self.val, pos-idx_of_pos)
        return QaryString(self.q, val)
    
    def _substitute(self, pos, symbol, idx_of_pos=0):
        val = np.copy(self.val)
        val[pos-idx_of_pos] = symbol
        return QaryString(self.q, val)
    
    @property
    def syndrome(self):
        return syndrome(self.val)
    
    @property
    def signature(self):
        return QaryString(q=2, val=signature(self.val))
    
    @property
    def parity_check(self):
        assert self.q == 2
        return parity_check(self.val)
    
    @property
    def sum(self):
        return np.sum(self.val)
    
    @property
    def marker(self):
        r = 1 if self.val[-1] == 0 else 0
        r = QaryString(self.q, val=np.array([r]))
        return r.concatenate([r])
    
    """
    Convenience functions for encoding as k-sum-balanced strings
    """
    @property
    def is_sum_balanced(self):
        return is_sum_balanced(self.val, self.q)
    
    @property
    def as_binary_matrix(self):
        logq = np.ceil(np.log(self.q) / np.log(2)).astype(self.val.dtype)
        m = np.zeros([self.length, logq], dtype=self.val.dtype)
        for i in range(self.length):
            binary = list(bin(self.val[i])[2:])
            binary = [0] * (logq - len(binary)) + binary            
            m[i,:] = np.array(binary)
        return self.q, m
    
    @staticmethod
    def from_binary_matrix(q, m):
        val = [QaryString(2, m[i]).asint() for i in range(m.shape[0])]
        return QaryString(q, val)