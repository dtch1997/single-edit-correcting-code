# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 00:55:22 2020

@author: Daniel Tan
"""

import numpy as np
import util

class QaryString:
    
    def __init__(self, q : int = 4, 
                       val = np.zeros(shape=0, dtype=np.int8)):
        self.q = q
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
            sum += self.val[i]
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
        return np.ceil(np.log(n) / np.log(self.q)).astype(int)
    
    def split(self, lengths):
        assert np.sum(np.array(lengths)) == self.length
        parts = []
        idx = 0
        for i, l in enumerate(lengths):
            parts.append(self[idx: idx + l])
            idx = idx + l
        return tuple(parts)
        
    def __len__(self):
        return self.val.shape[0]
    
    def __str__(self):
        return f"q = {self.q}, val = {self.val}"
    
    def __getitem__(self, key):
        return QaryString(self.q, self.val.__getitem__(key))
    
    def __setitem__(self, key, item):
        self.val.__setitem__(key, item)
        
    def __delitem__(self, item):
        self.val.__delitem__(item)
        
    def __eq__(self, other):
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
            val = np.delete(val, pos)
            symbol = None
        return QaryString(self.q, val), mtype, pos, symbol
         
    
    @property
    def syndrome(self):
        return util.syndrome(self.val)
    
    @property
    def signature(self):
        return QaryString(q=2, val=util.signature(self.val))
    
    @property
    def parity_check(self):
        assert self.q == 2
        return util.parity_check(self.val)
    
    @property
    def sum(self):
        return np.sum(self.val)
    
    @property
    def marker(self):
        r = 1 if self.val[-1] == 0 else 0
        r = QaryString(self.q, val=np.array([r]))
        return r.concatenate([r])
    
if __name__ == "__main__":
    for i in range(100):
        x = QaryString(val = np.zeros(2))
        xm, mtype = x.mutate()
        print(xm, mtype)
    