# Single Edit Correcting Code

Implements an order-optimal code for correcting single edits in a quartenary alphabet, described [here](https://arxiv.org/pdf/1910.06501.pdf)

## Encoding Scheme

The encoder needs to encode an arbitrary string **s** to a k-sum-balanced string **x**, then appends various checksums to **x** to create a final string **x_enc**, which is sent over some noisy channel.  

## Status


### Encoder
- [x] Implement utility code for working with q-ary strings. 
- [x] Implement code for encoding **x** to **x_enc**
- [x] Implement code for decoding a single substitution in **x_enc**
- [x] Implement code for decoding a single edit in **x_enc**
- [x] Implement code for decoding a single insertion in **x_enc**
- [x] Implement code for encoding **s** to **x**
- [x] Implement code for decoding **x** to **s** 

### Miscellaneous
- [ ] Update README with documentation on how to use code
- [ ] Update README with explanation of the encoding scheme
