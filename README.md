# Single Edit Correcting Code

An implementation of an order-optimal code for correcting single edits in a quartenary alphabet, first described [in this paper](https://arxiv.org/pdf/1910.06501.pdf)

## Getting Started

First, install requirements. We only require Numpy and Scipy. The versions suggested are not hard requirements; they are only a suggestion of a working version.
```
pip install -r requirements.txt
```
See `example.py` for sample usage of the code. 

## Encoding Scheme

The encoder needs to encode an arbitrary string **s** to a k-sum-balanced string **x**, then appends various checksums to **x** to create a final string **x_enc**, which is sent over some noisy channel. 

The main innovation is the k-sum-balanced string, which is a core component of the code construction. As discussed in Lemma 29 of [the paper](https://arxiv.org/pdf/1910.06501.pdf), it enables easy localization of insertion and deletion errors. This localization means that we can rely on an existing code known as a [shifted-VT code](https://arxiv.org/pdf/1602.06820.pdf) to find the insertion or deletion. 

## Contact

Questions about the code's theoretical guarantees can be directed to hmkiah@ntu.edu.sg  
Questions about bugs or implementation errors can be directed to dtch1997@stanford.edu
