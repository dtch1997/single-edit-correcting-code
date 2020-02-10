# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 23:28:57 2020

@author: Daniel Tan
"""

import numpy as np
import util
from qary_string import QaryString
from typing import List

class SingleEditCode:
    def __init__(self):
        pass
    
    def encode(self, x):
        """
        Parameters
        ----------
        x : balanced QaryString of length (n)

        Returns
        -------
        x_enc : QaryString of length (N)

        """
        x_enc = None
        n = x.length
        
        k = 72 * np.ceil(np.log(n) / np.log(2))
        P = 20*k
        
        a = x.syndrome % (4*n + 1)
        b = x.signature.syndrome % P
        c = x.signature.parity_check 
        d = x.sum % 7
        M = x.marker
        
        R1 = x.fromint(a).pad_to(x.bitlen(4*n+1))
        R2 = x.fromint(b).pad_to(x.bitlen(P))
        R3 = x.fromint(c).pad_to(1)
        R4 = x.fromint(d).pad_to(2)
        
        x_enc = x.concatenate([M, R1, R2, R3, R4])
        return x_enc, n, x_enc.length
    
    def decode(self, x_enc, n, N, verbose=False):
        """
        n:  int. Length of un-encoded message. 
        N:  int. Length of un-mutated encoded message. 
            
        x_enc: Encoded message, possibly mutated. 
        """
        assert N-1 <= x_enc.length <= N+1
        if x_enc.length == N:
            return self._decode_substitution(x_enc, n, verbose)
        elif x_enc.length == N+1:
            return self._decode_insertion(x_enc, n)
        elif x_enc.length == N-1:
            return self._decode_deletion(x_enc, n)
        # Should never reach this part
        raise Exception("x_enc has invalid length in decode()")
    
    def _decode_substitution(self, x_enc, n, verbose):
        k = 72 * np.ceil(np.log(n) / np.log(2))
        P = 20*k
        lengths = [n, 2, x_enc.bitlen(4*n+1), x_enc.bitlen(P), 1, 2]
        
        xp, Mp, R1p, R2p, R3p, R4p = x_enc.split(lengths)
        if Mp[0] != Mp[1]:
            return xp

        ap = xp.syndrome % (4*n + 1)
        dp = xp.sum % 7
        
        if R4p.asint() == dp:
            return xp
        if R1p.asint() == ap:
            return xp
        val_change = dp - R4p.asint()
        if val_change <= -x_enc.q: val_change += 7
        if val_change >= x_enc.q: val_change -= 7
        sig_change = ap - R1p.asint()
        for j in range(1, n+1):
            if (sig_change - j * val_change) % (4*n+1) == 0:
                break
            if j == n:
                raise Exception("j not found")
       
        if verbose:
            print(f"Detected a substitution at array index {j-1}")
            print(f"Increase in value: {val_change}")
        
        x = QaryString(xp.q, xp.val)
        x.val[j-1] = x.val[j-1] - val_change
        return x
    
    def _decode_deletion(self, x_enc):
        pass
    
    def _decode_insertion(self, x_enc):
        pass
    
if __name__ == "__main__":
    x = QaryString(4, np.array([1,1,2,2,1,1,2,2]))
    code = SingleEditCode()
    x_enc, n, N = code.encode(x)
    x_enc_p = QaryString(x_enc.q, x_enc.val)
    x_enc_p[2] = 1
    x_pred = code.decode(x_enc_p, n, N)
    print(x_pred)
    print(x)