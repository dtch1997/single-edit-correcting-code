# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 23:52:06 2020

@author: Daniel Tan
"""


from dnacode import *
from find_sum_balanced import *
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


        
if __name__ == "__main__":
    unittest.main()
