# Single Edit Correcting Code

Implements an order-optimal code for correcting single edits in a quartenary alphabet, described [here](https://arxiv.org/pdf/1910.06501.pdf)

## Encoding Scheme

The encoder needs to encode an arbitrary string **x** to a k-sum-balanced string **y**, then appends various checksums to **y** to create a final string **z**, which is sent over some noisy channel.  

## Status

- [x] Implement utility code for working with q-ary strings. 
- [x] Implement code for encoding **y** to **z**
- [x] Implement code for decoding a single substitution in **z**
- [ ] Implement code for decoding a single edit in **z**
- [ ] Implement code for decoding a single insertion in **z**
- [ ] Implement code for encoding **x** to **y**
- [ ] Implement code for encoding **y** to **z** 
