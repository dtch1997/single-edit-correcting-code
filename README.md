# Single Edit Correcting Code

An implementation of an order-optimal code for correcting single edits (insertion, deletion, or substitution) in a quartenary alphabet, first described in section IV-B of the paper cited below. We thank the authors of this work for their help with the implementation.

Cai, Kui, Yeow Meng Chee, Ryan Gabrys, Han Mao Kiah, and Tuan Thanh Nguyen. "Optimal codes correcting a single indel/edit for DNA-based data storage." _arXiv preprint arXiv:1910.06501_ (2019). [arXiv](https://arxiv.org/pdf/1910.06501.pdf).

This work was done by Daniel Tan as a course project for EE 276 (Information Theory) - Winter 2020 at Stanford University, taught by Prof. Tsachy Weissman, and this project was mentored by TA Shubham Chandak. The brief report for the project describing the links to DNA storage is available [here](https://theinformaticists.com/2020/03/25/an-order-optimal-edit-correcting-code/).

See the repository [shubhamchandak94/VT_codes](https://github.com/shubhamchandak94/VT_codes) for the implementation of asymptotically optimal Varshamov-Tenengolts (VT) codes that can correct a single insertion/deletion for general alphabets, and optionally also correct a single substitution for the binary case.

## Getting Started

Clone the repository:
```
git clone git@github.com:dtch1997/single-edit-correcting-code.git
cd single-edit-correcting-code
```
Install requirements. We require Python3, `numpy` and `scipy`, as well as `nose` for testing.
```
pip install -r requirements.txt
```
Run the test suite.
```
nosetests
```
See `example.py` for sample usage of the code with detailed comments. Run as `python example.py` for a simple demo.

The organization of the library:
- [single_edit_code/single_edit_code.py](single_edit_code/single_edit_code.py): the main library with code creation, encoding and decoding.
- [single_edit_code/qary_string.py](single_edit_code/qary_string.py): class for manipulating q-ary strings.
- [single_edit_code/svt_code.py](single_edit_code/svt_code.py): shifted-VT code used internally for correcting localized errors.
- [single_edit_code/sum_balanced_code.py](single_edit_code/sum_balanced_code.py): sum-balanced codes used internally for the encoding/decoding.
- [single_edit_code/util.py](single_edit_code/util.py): Utility functions.

The implementation currently has certain limitations (not necessarily shared by the paper):
- message to be encoded must already be quaternary string, binary strings are not supported.
- Certain encoding/decoding related parameters (such as length of intermediate sum-balanced string and the length of the codeword) can be computed only by performing the encoding. This is inconvenient because computing the optimal message length for a given codeword length is non-trivial, and these additional parameters must be provided to the decoder. In addition, some of the parameters slightly depend on the specific message, and hence need to be stored in addition to the codeword in the current implementation.


## Encoding Scheme

In the following section we provide a high-level overview of how the encoding and decoding process works, as implemented in this repository. Please refer to the original paper for a more detailed and rigorous discussion.

The encoder encodes a string **s** to a k-sum-balanced string **x** (Definition 24), then appends various checksums to **x** to create a final string **x_enc**, which is sent over some noisy channel. The decoder receives a possibly-corrupted message **x_enc_m** and recovers the original string **s**. All strings are assumed to be quaternary - i.e. over an alphabet of four symbols 0,1,2,3.

The main innovation of this paper is the k-sum-balanced string, which is a core component of the code construction. As discussed in Lemma 29 of [the paper](https://arxiv.org/pdf/1910.06501.pdf), it enables easy localization of insertion and deletion errors. This localization means that we can rely on an existing code known as a [shifted-VT code](https://arxiv.org/pdf/1602.06820.pdf) to find the insertion or deletion. Substitution errors are easier to detect and do not require this localization.

### Creating a k-sum-balanced string

The original paper proves that an arbitrary string **s** can be encoded into a k-sum-balanced string **x** that is only one symbol longer for a sufficiently large choice of **k**. The following high-level procedure is used:

1. Move a sliding window of length **k** over **s**. At each position, we look for "forbidden words" which are not k-sum-balanced.
2. If a **k**-slice **w** is not k-sum-balanced, we encode it with a tuple of `rep_w = (ind_w, pos, 3)` where:
    * **ind_w** is the index of **w** in an enumeration of all forbidden words **w**
    * **pos** describes the position where **w** was found
    * 3 is an end-of-block marker symbol.
3. We can then remove **w** from the sequence and append `rep_w` to the tail of the sequence.

This process can be repeated until no forbidden words remain in the sequence.  

The decoding proceeds in reverse: We read in blocks of `(ind_w, pos, 3)` from the tail of the sequence, recover the original **w** and insert it at the appropriate location.

### Computing index of a forbidden word

The naive method of enumerating all possible forbidden words is an `O(4^k)` operation which is unacceptably long for all reasonable values of **k**. We use an index built on the [combinatorial number system](https://en.wikipedia.org/wiki/Combinatorial_number_system) to compute indices in `O(k^2)` time given `O(k^3)` precomputation time. This is a good choice because there are exponentially many possible forbidden words but only a few of them will ever be encountered in a typical message.

## Contact

Questions about the code's theoretical guarantees can be directed to Han Mao Kiah (author of previously cited paper).
Questions about bugs or implementation errors can be raised as GitHub issues.
