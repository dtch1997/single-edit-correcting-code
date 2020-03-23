# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 23:31:11 2020

@author: Daniel Tan
"""

import numpy as np
import single_edit_code as sec

x = sec.QaryString(4, np.random.randint(0, high=2, size=40))
code = sec.SingleEditCode(16)
x_enc, n, N, l = code.encode(x)
x_enc_m, _, pos, symbol = x_enc.mutate(mtype="insert")
x_pred = code.decode(x_enc_m, n, N, l, verbose=False)

assert x_pred == x
            