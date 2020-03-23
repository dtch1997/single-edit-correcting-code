# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 23:39:55 2020

@author: Daniel Tan
"""
import numpy as np
import single_edit_code as sec
import pdb

k = 36

code = sec.SumBalancedCode(k)
for i in range(1000):
    s = sec.QaryString(q=4, val=np.random.randint(low=0, high=4, size=100))
    
    try: 
        x, l = code.encode(s)
    except Exception as e:
        print(e)
        print(code.sumpair2bucket)
        raise Exception()
        
    try:
        assert sec.is_k_sum_balanced(x.val, k)
    except Exception as e:
        print("Encoded string was not sum balanced")
        print(s)
        print(x)
        raise Exception()
        
    s_pred = code.decode(x, l)
    try: 
        assert s == s_pred
    except Exception as e:
        print("Decoded string was wrong")
        print(s)
        print(s_pred)
        print(code.k)
        print(code.bucket2startidx)
        raise Exception()
        
        
    print(f"Success {i}")
    