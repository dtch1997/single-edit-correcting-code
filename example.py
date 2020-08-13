# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 23:31:11 2020

@author: Daniel Tan
"""

import numpy as np
import single_edit_code as sec

# generate random quaternary string of length 40
x = sec.QaryString(4, np.random.randint(0, high=3, size=160))

# print quaternary string to visualize
print(x)

# the value of k determines the redundancy of the code,
# it should be roughly 36*log_4(codeword_length)
# Since codeword_length might not be known a priori, we can set to
# roughly 36*log_4(message_length) which is 36*log_4(160) ~ 132 in our case.
# It is possible to use lower values of k, but then the redundancy is higher. Note
# that higher k = more precomputation and slower encoding/decoding.

# this command constructs the code, this does some precomputation such
# as building index for forbidden words, that reduces the encoding time later.
code = sec.SingleEditCode(k=132)

# perform the encoding, returning the encoded quaternary string string + parameters needed during decoding.
x_enc, n, N, l = code.encode(x)
# print quaternary string to visualize
print('x_enc', x_enc)
print('length of message: ', l)
print('length of codeword: ', N)
print('length of intermediate sum-balanced codeword: ', n)

# now we mutate the codeword in various ways and verify that decoding succeeds

# test substitution
print('**Test substitution error**')
x_enc_m, _, pos, symbol = x_enc.mutate(mtype="substitute")
print('Substitution from', x_enc.val[pos], 'to', x_enc_m.val[pos], 'at position', pos)
assert len(x_enc_m) == N
# note that decoder requires the values l, N and n, so these need to be retained in
# addition to the codeword itself
x_pred = code.decode(x_enc_m, n, N, l, verbose=False)
assert x_pred == x
print("Decoding succeeded")

print('**Test insertion error**')
x_enc_m, _, pos, symbol = x_enc.mutate(mtype="insert")
print('Insertion of', symbol, 'at position', pos)
assert len(x_enc_m) == N+1
# note that decoder requires the values l, N and n, so these need to be retained in
# addition to the codeword itself
x_pred = code.decode(x_enc_m, n, N, l, verbose=False)
assert x_pred == x
print("Decoding succeeded")

print('**Test deletion error**')
x_enc_m, _, pos, symbol = x_enc.mutate(mtype="delete")
print('Deletion at position', pos)
assert len(x_enc_m) == N-1
# note that decoder requires the values l, N and n, so these need to be retained in
# addition to the codeword itself
x_pred = code.decode(x_enc_m, n, N, l, verbose=False)
assert x_pred == x
print("Decoding succeeded")
