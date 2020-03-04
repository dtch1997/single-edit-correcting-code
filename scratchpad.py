# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 22:08:18 2020

@author: Daniel Tan
"""


import numpy as np
from util import *
from qary_string import QaryString
from itertools import product

k=4

count = 0
for seq in product([0,1,2,3], repeat=k):
    x = np.array(seq)
    if not is_k_sum_balanced(x, k):
        count += 1
        print(x)
print(f"Number of forbidden words: {count} / {4 ** k}")
        