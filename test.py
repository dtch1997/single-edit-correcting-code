# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 23:52:06 2020

@author: Daniel Tan
"""


from dnacode import *
from util import *
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
        
class TestEncode(unittest.TestCase):
    def test_seq12(self):
        """
        Run tests on sequences of 1s and 2s, which will (almost) always be sum-balanced. 
        """
        for i in range(1000):
            x = QaryString(4, np.random.randint(1, 3, 256))
            code = SingleEditCode()
            x_enc, n, N = code.encode(x)
            x_enc_m, mtype, pos, symbol = x_enc.mutate(mtype="substitute")
            x_pred = code.decode(x_enc_m, n, N)
            self.assertTrue(x_pred == x)
            
    def test1_fixed(self):
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
        
        
if __name__ == "__main__":
    unittest.main()
