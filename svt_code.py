# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 22:19:12 2020

@author: Daniel Tan
"""


from qary_string import QaryString

def first_index_k_zeros_left(qstr, k, P):
    """
    For a binary string qstr, return the first index of q with k (mod P) zeros to the left.
    Return: index in [0, qstr.length]
    """
    num_zeros_left = 0
    for j in range(qstr.length+1):
        if (num_zeros_left - k) % P == 0:
            return j
        if j == qstr.length: 
            raise Exception("Should never reach here")
        if qstr[j] == 0:
            num_zeros_left += 1
            
def first_index_k_ones_right(qstr, k, P):
    """
    For a binary string qstr, return the first index of q with k (mod P) ones to the right.
    Return: index in [0, qstr.length]
    """
    num_ones_right = 0
    for j in range(qstr.length, -1, -1):
        if (num_ones_right - k) % P == 0:
            return j
        if j == 0: 
            raise Exception("Should never reach here")
        if qstr[j-1] == 1:
            num_ones_right += 1

class SVTCode:
    """
    Implemented based on Algorithm 2 in 
    https://arxiv.org/pdf/1602.06820.pdf
    """
        
    def __init__(self):
        pass
    
    def encode(self, x):
        """
        x: QaryString with q=2, length=n
        """
        x_enc = None
        return x_enc
    
    def decode_deletion(self, y, a, u, P, delval, verbose=False):
        """
        y: A QaryString, q=2 with one symbol deleted.  
        a: VT syndrome of the original string. 
        u: first possible position of the deletion, 1-indexed
        P: upper bound on number of positions where deletion could have occurred. 
        delval: deleted symbol. One of 0,1
        """
        P = min(P, y.length - u + 2)
        lengths = [u-1, P-1, y.length-u-P+2]
        if verbose: print(lengths)
        yhead, yhat, ytail = y.split(lengths)
        if verbose: print(f"yhead = {yhead} \nyhat = {yhat} \nytail = {ytail}")        
        if verbose: print(f"Running SVT decode with u = {u}, P={P}, delval={delval}, a={a}")
        
        ap = (y.syndrome + ytail.sum) % P
        if verbose: print(f"Augmented weighted sum of {ap}")
        delta = (a - ap) % P
        
        if verbose: print(f"delta = {delta}")
        
        if delval == 0:
            delpos = first_index_k_ones_right(yhat, delta, P)
            # Zero indexed
            # First position to the left of delta ones, not necessarily consecutive
            """
            num_ones_right = 0
            for k in range(yhat.length, -1, -1):
                if num_ones_right == delta:
                    delpos = k
                    break
                if yhat[k-1] == 1:
                    num_ones_right += 1
            if delpos is None: raise Exception("delval = 0 should not reach here")
            """
        if delval == 1:
            num_zeros_left = (delta - u - yhat.sum) % P
            if verbose: print(f"finding first position to the right of " 
                              f"{delta} - {u} - {yhat.sum} (mod {P}) = "                           
                              f"{num_zeros_left} zeros")
            delpos = first_index_k_zeros_left(yhat, delta - u - yhat.sum, P)
            # Zero indexed
            # First position to the right of delta - mu - wt(yhat) zeros
            
        if verbose: print(f"delpos = {delpos}")
        yhat = yhat._insert(delpos, delval, idx_of_pos = 0)
        if verbose: print(f"New yhat: {yhat}")
        return yhead.concatenate(yhat).concatenate(ytail)
    
    def decode_insertion(self, y, a, u, P, insval, verbose=False):
        """
        y: A QaryString, q=2 with one symbol inserted. 
        a: VT syndrome of the original string. 
        u: first possible position of the insertion, 1-indexed
        P: upper bound on number of positions where insertion could have occurred. 
        insval: inserted symbol. One of 0,1
        """
        P = min(P, y.length - u + 1)
        lengths = [u-1, P, y.length-u-P+1]
        if verbose: print(f"Splitting y according to lengths {lengths}")
        yhead, yhat, ytail = y.split(lengths)
        if verbose: print(f"yhead = {yhead} \nyhat = {yhat} \nytail = {ytail}")        
        if verbose: print(f"Running SVT decode with u = {u}, P={P}, insval={insval}, a={a}")
        
        ap = (y.syndrome - ytail.sum) % P
        if verbose: print(f"Augmented weighted sum of {ap}")
        delta = (ap - a) % P
        
        if verbose: print(f"delta = {delta}")
        
        if insval == 0:
            # A zero was inserted to the left of delta ones. 
            if verbose: print(f"Finding first position to the left of {delta} ones")
            inspos = first_index_k_ones_right(yhat, delta, P)
            # Minus one because the above function was written for the deletion case. 
            inspos = inspos - 1
        if insval == 1:
            # A one was inserted to the right of (delta - u - wy(yhat)) zeros. 
            if verbose: print(f"Finding first position to the right of {delta} - "
                  f"{u} - {yhat.sum} + 1 (mod {P}) = "
                  f"{(delta - u - yhat.sum + 1) % P} zeros")
            inspos = first_index_k_zeros_left(yhat, delta - u - yhat.sum + 1, P)

        yhat = yhat._delete(inspos)
        return yhead.concatenate(yhat).concatenate(ytail)
    
    
    
if __name__ == "__main__":
    y = QaryString(q=2, val=[0,1]*5)
    pos = 4
    symbol = 1
    y_ins = y._insert(pos, symbol, idx_of_pos=1)
    print(f"Inserted symbol {symbol} at {pos}")
    print(f"Mutated sequence {y_ins}")
    y_pred = SVTCode().decode_insertion(y_ins, y.syndrome, pos, 3, symbol, verbose=True)
    print(f"Predicted {y_pred}")
        