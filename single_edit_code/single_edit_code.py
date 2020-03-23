# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 23:28:57 2020

@author: Daniel Tan
"""

import numpy as np
from .qary_string import QaryString
from .svt_code import SVTCode
from .sum_balanced_code import SumBalancedCode
from typing import List

class SingleEditCode:
    def __init__(self, k: int = 64):
        if type(k) is not int:
            raise Exception("SingleEditCode requires integer k")
        self.k = k
        self.sbcode = SumBalancedCode(k)
    
    def encode(self, x):
        """
        Parameters
        ----------
        x : QaryString of length (n)

        Returns
        -------
        x_enc : QaryString of length (N)

        """
        
        x, l = self.sbcode.encode(x)        
        x_enc = None
        n = x.length
        
        k = self.k
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
        return x_enc, n, x_enc.length, l
    
    def decode(self, x_enc, n, N, l, verbose=False):
        """
        n:  int. Length of un-encoded message. 
        N:  int. Length of un-mutated encoded message. 
            
        x_enc: Encoded message, possibly mutated. 
        """
        if verbose: print("Beginning decoding...")
        assert N-1 <= x_enc.length <= N+1
        if x_enc.length == N:
            x_dec_ksumbalanced = self._decode_substitution(x_enc, n, verbose)
        elif x_enc.length == N+1:
            x_dec_ksumbalanced = self._decode_insertion(x_enc, n, verbose)
        elif x_enc.length == N-1:
            x_dec_ksumbalanced = self._decode_deletion(x_enc, n, verbose)
        else: 
            raise Exception("x_enc has invalid length in decode()")
        return self.sbcode.decode(x_dec_ksumbalanced, l)
            
    
    def _decode_substitution(self, x_enc, n, verbose):
        k = self.k
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
    
    def _decode_deletion(self, x_enc, n, verbose):
        if verbose: 
            print("Deletion detected.")
            print("Computing parameters...")
        P = self._get_P(n)
        lengths = [n-1, 2, x_enc.bitlen(4*n+1), x_enc.bitlen(P), 1, 2]
        xp, Mp, R1p, R2p, R3p, R4p = x_enc.split(lengths)
        if Mp[0] != Mp[1]:
            return xp.concatenate(Mp[0])
        
        ap = xp.syndrome % (4*n + 1)
        cp = xp.signature.parity_check
        dp = xp.sum % 7
        
        xp_deleted_symbol = (R4p.asint() - dp) % 7
        if verbose: print(f"Deletion of symbol {xp_deleted_symbol} detected.")
        
        Js = []
        # 1-indexed positions where the deletion could have happened. 
        for j in range(1, n+1):
            if (ap + j*xp_deleted_symbol + xp[j-1:].sum) % (4*n+1) == R1p.asint() % (4*n+1):
                Js.append(j)
        if verbose: print(f"Possible locations: {Js}")
        
        Jp = set()
        for j in Js:
            Jj = set()
            xp_j = xp._insert(j, xp_deleted_symbol, idx_of_pos=1)
            for jp in range(1,xp_j.signature.length +1):
                if xp_j.signature._delete(jp, idx_of_pos=1) == xp.signature:
                    Jj.add(jp)
            Jp = Jp.union(Jj)
            
        if verbose: print(f"Possible deletion locations: {Jp}")
        u = min(Jp)
        if verbose: print(f"P = {P}")
        sig_deleted_symbol = (cp - R3p.asint()) % 2
        sig = SVTCode().decode_deletion(xp.signature, R2p.asint(), u, P, sig_deleted_symbol, verbose=verbose)
        
        for j in Js:
            x_candidate = xp._insert(j, xp_deleted_symbol, idx_of_pos=1)
            if x_candidate.signature == sig:
                return x_candidate
            
        raise Exception("SingleEditCode could not decode the given string")
        
                    
    def _decode_insertion(self, x_enc, n, verbose):
        if verbose: 
            print("Insertion detected.")
            print("Computing parameters...")
        P = self._get_P(n)
        lengths = [n+1, 2, x_enc.bitlen(4*n+1), x_enc.bitlen(P), 1, 2]
        xp, Mp, R1p, R2p, R3p, R4p = x_enc.split(lengths)
        if xp[-1] == Mp[0]:
            return xp[:-1]
        
        ap = xp.syndrome % (4*n + 1)
        cp = xp.signature.parity_check
        dp = xp.sum % 7
        xp_inserted_symbol = (dp - R4p.asint()) % 7
        if verbose: print(f"Insertion of symbol {xp_inserted_symbol} detected.")
        
        Js = []
        # 1-indexed positions where insertion could have occured
        for j in range(1, n+2):
            if xp.val[j-1] == xp_inserted_symbol and \
            (ap - j*xp_inserted_symbol - xp[j:].sum) % (4*n+1) == R1p.asint() % (4*n+1):
                Js.append(j)
        if verbose: print(f"Possible locations: {Js}")
        
        Jp = set()
        for j in Js:
            Jj = set()
            xp_j = xp._delete(j, idx_of_pos=1)
            for jp in range(1,xp.signature.length +1):
                if xp_j.signature == xp.signature._delete(jp, idx_of_pos=1):
                    Jj.add(jp)
            Jp = Jp.union(Jj)
        if verbose: print(f"Possible insertion locations: {Jp}")
        u = min(Jp)
        if verbose: print(f"P = {P}")
        
        sig_inserted_symbol = (cp - R3p.asint()) % 2
        sig = SVTCode().decode_insertion(xp.signature, R2p.asint(), u, P, sig_inserted_symbol, verbose=verbose)
        
        for j in Js:
            x_candidate = xp._delete(j, idx_of_pos=1)
            if x_candidate.signature == sig:
                return x_candidate
            
        raise Exception("SingleEditCode could not decode the given string")

    def _get_k(self, n):
        return self.k
        # return 72 * np.ceil(np.log(n) / np.log(2))
    
    def _get_P(self, n):
        return 20*self._get_k(n)
            