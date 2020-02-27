# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 23:52:06 2020

@author: Daniel Tan
"""


from single_edit_code import *
from util import *
from svt_code import *
import unittest
import numpy as np

class TestSyndrome(unittest.TestCase):
    def test1(self):
        x = np.array([0,0,1])
        self.assertEqual(syndrome(x), 3)
        
class TestSignature(unittest.TestCase):
    def test1(self):
        x = np.array([3,2,5])
        result = np.all(signature(x) == np.array([False, True]))
        self.assertTrue(result)
        
class TestSumBalanced(unittest.TestCase):
    def test1_sum_balanced(self):
        x = np.array([0,0,0,1,1,1])
        self.assertFalse(is_sum_balanced(x))
        self.assertTrue(is_sum_balanced(np.arange(4)))

class TestkSumBalanced(unittest.TestCase):
    def test1_k_sum_balanced(self):
        self.assertFalse(is_k_sum_balanced(np.arange(4), 1))
        self.assertTrue(is_k_sum_balanced(np.array([1,1,2,2,1,1,2,2]), 3))
        
class TestSingleEditCode(unittest.TestCase):
    def test_substitute_seq12(self):
        """
        Run tests on sequences of 1s and 2s, which will (almost) always be sum-balanced. 
        """
        for i in range(10):
            x = QaryString(4, np.random.randint(1, 3, 2048))
            code = SingleEditCode()
            x_enc, n, N = code.encode(x)
            x_enc_m, mtype, pos, symbol = x_enc.mutate(mtype="substitute")
            x_pred = code.decode(x_enc_m, n, N)
            self.assertTrue(x_pred == x)
            
    def test1_substitute_fixed(self):
        """
        Run tests on a fixed sequence which we know to be sum balanced
        """
        for i in range(1000):
            x = QaryString(4, np.array([0,1,2,3]*10))
            code = SingleEditCode()
            x_enc, n, N = code.encode(x)
            x_enc_m, mtype, pos, symbol = x_enc.mutate(mtype="substitute")
            x_pred = code.decode(x_enc_m, n, N)
            self.assertTrue(x_pred == x)
             
    def test1_delete_fixed(self):
        # Run tests on a fixed sequence which we know to be sum balanced
        for i in range(1000):
            x = QaryString(4, np.array([0,1,2,3]*10))
            code = SingleEditCode()
            x_enc, n, N = code.encode(x)
            x_enc_m, mtype, pos, symbol = x_enc.mutate(mtype="delete")
            x_pred = code.decode(x_enc_m, n, N)
            self.assertTrue(x_pred == x)

    def test_delete_seq12(self):
        # Run tests on sequences of 1s and 2s, which will (almost) always be sum-balanced. 
        for i in range(10):
            x = QaryString(4, np.random.randint(1, 3, 2048))
            code = SingleEditCode()
            x_enc, n, N = code.encode(x)
            x_enc_m, mtype, pos, symbol = x_enc.mutate(mtype="delete")
            x_pred = code.decode(x_enc_m, n, N)
            self.assertTrue(x_pred == x)
            
    def test1_insert_fixed(self):
        # Run tests on a fixed sequence which we know to be sum balanced
        for i in range(1000):
            x = QaryString(4, np.array([0,1,2,3]*10))
            code = SingleEditCode()
            x_enc, n, N = code.encode(x)
            x_enc_m, mtype, pos, symbol = x_enc.mutate(mtype="insert")
            x_pred = code.decode(x_enc_m, n, N)
            self.assertTrue(x_pred == x)

    def test_insert_seq12(self):
        # Run tests on sequences of 1s and 2s, which will (almost) always be sum-balanced. 
        for i in range(10):
            x = QaryString(4, np.random.randint(1, 3, 2048))
            code = SingleEditCode()
            x_enc, n, N = code.encode(x)
            x_enc_m, mtype, pos, symbol = x_enc.mutate(mtype="insert")
            x_pred = code.decode(x_enc_m, n, N)
            self.assertTrue(x_pred == x)
            
    def test_comprehensive(self):
        for i in range(1000):
            x = QaryString(4, np.random.randint(4, 1024))
            code = SingleEditCode()
            x_enc, n, N = code.encode(x)
            x_enc_m, mtype, pos, symbol = x_enc.mutate(mtype="insert")
            x_pred = code.decode(x_enc_m, n, N)
            self.assertTrue(x_pred == x)
            

class TestSVTCode(unittest.TestCase):
    def test1_delete_fixed_short(self):
        y = QaryString(q=2, val=[0,1,0,1,0,1,0,1])
        symbol = 0
        pos = 4
        yp = y._delete(pos, idx_of_pos = 0)
        code = SVTCode()
        y_pred = code.decode_deletion(yp, y.syndrome, pos+1, 3, symbol, verbose=False)
        self.assertEqual(y_pred, y)
        
        
    def test1_delete_fixed_long(self):
        y = QaryString(q=2, val=[0,1]*1000)
        symbol = 0
        pos = 4
        yp = y._delete(pos, idx_of_pos = 0)
        code = SVTCode()
        y_pred = code.decode_deletion(yp, y.syndrome, pos+1, 3, symbol, verbose=False)
        self.assertEqual(y_pred, y)
        

    def test1_delete_random_short(self):
        for i in range(1000):
            y = QaryString(2, np.random.randint(0, 2, 10))
            code = SVTCode()
            yp, _, pos, symbol = y.mutate(mtype="delete")
            u = np.random.randint(low=1, high=pos+2)
            y_pred = code.decode_deletion(yp, y.syndrome, u, 10, symbol, verbose=False)
            self.assertTrue(y_pred == y)
            
    def test1_insert_fixed_short(self):
        y = QaryString(q=2, val=[0,1,0,1,0,1,0,1])
        symbol = 0
        pos = 4
        yp = y._insert(pos, symbol, idx_of_pos = 0)
        code = SVTCode()
        y_pred = code.decode_insertion(yp, y.syndrome, pos+1, 3, symbol, verbose=False)
        self.assertEqual(y_pred, y)
        
        
    def test1_insert_fixed_long(self):
        y = QaryString(q=2, val=[0,1]*1000)
        symbol = 0
        pos = 4
        yp = y._insert(pos, symbol, idx_of_pos = 0)
        code = SVTCode()
        y_pred = code.decode_insertion(yp, y.syndrome, pos+1, 3, symbol, verbose=False)
        self.assertEqual(y_pred, y)
        

    def test1_insert_random_short(self):
        for i in range(1000):
            y = QaryString(2, np.random.randint(0, 2, 10))
            code = SVTCode()
            yp, _, pos, symbol = y.mutate(mtype="insert")
            u = np.random.randint(low=1, high=pos+2)
            y_pred = code.decode_insertion(yp, y.syndrome, u, 10, symbol, verbose=False)
            self.assertTrue(y_pred == y)
            
if __name__ == "__main__":
    unittest.main()
