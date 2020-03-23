# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 23:52:06 2020

@author: Daniel Tan
"""

import single_edit_code as sec
import unittest
import numpy as np

class TestSyndrome(unittest.TestCase):
    def test1(self):
        x = np.array([0,0,1])
        self.assertEqual(sec.util.syndrome(x), 3)
        
class TestSignature(unittest.TestCase):
    def test1(self):
        x = np.array([3,2,5])
        result = np.all(sec.util.signature(x) == np.array([False, True]))
        self.assertTrue(result)
        
class TestSumBalanced(unittest.TestCase):
    def test1_sum_balanced(self):
        x = np.array([0,0,0,1,1,1])
        self.assertFalse(sec.util.is_sum_balanced(x))
        self.assertTrue(sec.util.is_sum_balanced(np.arange(4)))

class TestkSumBalanced(unittest.TestCase):
    def test1_k_sum_balanced(self):
        self.assertFalse(sec.util.is_k_sum_balanced(np.arange(4), 1))
        self.assertTrue(sec.util.is_k_sum_balanced(np.array([1,1,2,2,1,1,2,2]), 3))
        
class TestSingleEditCode(unittest.TestCase):            
    def test_comprehensive(self):
        code = sec.SingleEditCode()
        for i in range(1000):
            x = sec.QaryString(4, np.random.randint(4, 1024))
            x_enc, n, N, l = code.encode(x)
            x_enc_m, mtype, pos, symbol = x_enc.mutate(mtype="insert")
            x_pred = code.decode(x_enc_m, n, N, l)
            self.assertTrue(x_pred == x)
            
    def test_very_large_k(self):
        code = sec.SingleEditCode(512)
        for i in range(1000):
            x = sec.QaryString(4, np.random.randint(4, 1024))
            x_enc, n, N, l = code.encode(x)
            x_enc_m, mtype, pos, symbol = x_enc.mutate(mtype="insert")
            x_pred = code.decode(x_enc_m, n, N, l)
            self.assertTrue(x_pred == x)
            

class TestSVTCode(unittest.TestCase):
    def test1_delete_fixed_short(self):
        y = sec.QaryString(q=2, val=[0,1,0,1,0,1,0,1])
        symbol = 0
        pos = 4
        yp = y._delete(pos, idx_of_pos = 0)
        code = sec.SVTCode()
        y_pred = code.decode_deletion(yp, y.syndrome, pos+1, 3, symbol, verbose=False)
        self.assertEqual(y_pred, y)
        
        
    def test1_delete_fixed_long(self):
        y = sec.QaryString(q=2, val=[0,1]*1000)
        symbol = 0
        pos = 4
        yp = y._delete(pos, idx_of_pos = 0)
        code = sec.SVTCode()
        y_pred = code.decode_deletion(yp, y.syndrome, pos+1, 3, symbol, verbose=False)
        self.assertEqual(y_pred, y)
        

    def test1_delete_random_short(self):
        for i in range(1000):
            y = sec.QaryString(2, np.random.randint(0, 2, 10))
            code = sec.SVTCode()
            yp, _, pos, symbol = y.mutate(mtype="delete")
            u = np.random.randint(low=1, high=pos+2)
            y_pred = code.decode_deletion(yp, y.syndrome, u, 10, symbol, verbose=False)
            self.assertTrue(y_pred == y)
            
    def test1_insert_fixed_short(self):
        y = sec.QaryString(q=2, val=[0,1,0,1,0,1,0,1])
        symbol = 0
        pos = 4
        yp = y._insert(pos, symbol, idx_of_pos = 0)
        code = sec.SVTCode()
        y_pred = code.decode_insertion(yp, y.syndrome, pos+1, 3, symbol, verbose=False)
        self.assertEqual(y_pred, y)
        
        
    def test1_insert_fixed_long(self):
        y = sec.QaryString(q=2, val=[0,1]*1000)
        symbol = 0
        pos = 4
        yp = y._insert(pos, symbol, idx_of_pos = 0)
        code = sec.SVTCode()
        y_pred = code.decode_insertion(yp, y.syndrome, pos+1, 3, symbol, verbose=False)
        self.assertEqual(y_pred, y)
        

    def test1_insert_random_short(self):
        for i in range(1000):
            y = sec.QaryString(2, np.random.randint(0, 2, 10))
            code = sec.SVTCode()
            yp, _, pos, symbol = y.mutate(mtype="insert")
            u = np.random.randint(low=1, high=pos+2)
            y_pred = code.decode_insertion(yp, y.syndrome, u, 10, symbol, verbose=False)
            self.assertTrue(y_pred == y)
            
class TestCombinatorialBitstringEncoder(unittest.TestCase):
    def test1_random(self):          
        for i in range(1000):
            b = np.random.randint(low=0, high=2, size=64)
            k, w, index = sec.CombinatorialBitstringEncoder.encode(b)
            b_pred = sec.CombinatorialBitstringEncoder.decode(k,w,index)
            self.assertTrue(np.all(b_pred == b))
            
class TestSumBalancedCode(unittest.TestCase):
    k = 32
    code = sec.SumBalancedCode(k)
    for i in range(1000):
        s = sec.QaryString(q=4, val=np.random.randint(low=0, high=4, size=100))
        x, l = code.encode(s)
        assert sec.is_k_sum_balanced(x.val, k, q=4)
        s_pred = code.decode(x, l)
        assert s == s_pred
        
if __name__ == "__main__":
    unittest.main()
