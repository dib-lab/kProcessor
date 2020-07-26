# FAQ

<hr>

## **1. Hashing**

### 1.1 What are kDataFrame hashing modes?

> **Hashing modes are only applied in kDataFrameMQF and kDataFramePHMAP**

There are multiple hasing modes for the `kDataFrame` to use.
- 0: When setting the mode parameter to `0` that will define the hashing to be irreversible. By other words, once the kmer is hashed and inserted in the `kDataFrame` can't be reversed back to it's string representation.
- 1: Hashing mode `1` will define the hashing to be reversible, so that you can iterate over the `kDataFrame` kmers in their string representation. 

### 1.2 Why there are not user-selected hashing modes in `kDataFrameMAP`?

The only purpose of using the `kDataFrameMAP` is to store the kmers lexicographically. Consequently, no hashing is required, just store the kmers in their Two-bits representation.

### 1.3 What's the best hashing mode to use?

There's no best hashing mode, that depends totally on the application.

The irreversible mode is way faster than the reversible mode and requires less space. However, it's not as sensitive as the reversible mode in hashing the kmers because it stores the approximate hash value of the kmers with a user-controlled false-positive rate.
As the false-positive rate increases, the probability of hashing multiple kmers with the same hash value will increase.
Moreover, querying the kDataFrame with a similar but not stored kmer can show that the kmer already exist, but actually it's not. 

On the other hand, the reversible hashing mode is slower than the irreversible mode but allows the user to iterate over the kmers in their string representation.

---

## **2. Sequence Parsing Parameters**

In order to extract the kmers from a sequences file or string in kProcessor, you will need to pass a Python dictionary with the parsing parameters.

### 2.1 Parameters

Keys of the dictionary are *strings* describing the parameter type, values are integers for the parameters values

### 2.2 Modes

There are three modes for decoding the sequences substrings `kmers`, `skipmers`, `minimizers`

#### 2.2.1 kmers
    
- Description: Extracts the sequences substrings in the default popular mode "kmers".
- Parameters: 
    - "mode" : 1
    - k_size: total number of bases

#### 2.2.2 skipmers

- Description: A cyclic pattern of picked or skipped positions. [read more about skipmers](https://www.biorxiv.org/content/10.1101/179960v2)
- Parameters:
    - "mode" : 2
    - "k_size": total number of bases 
    - "m": used bases per cycle
    - "n": cycle length

#### 2.2.3 minimizers

- Description: short substrings that represents the sequence. [read more about minimizers](https://homolog.us/blogs/bioinfo/2017/10/25/intro-minimizer/)
- Parameters:
    - "mode" : 3
    - "k_size": total number of bases 
    - "w": window size

### 2.3 Examples

#### 2.3.1 Extract kmers with kmer size 31

```python
parse_params = {
    "mode" : 1,
    "k_size" : 31
}
```

#### 2.3.2 Extract skipmers with k = 10, m = 2, n = 3

```python
parse_params = {
    "mode" : 2,
    "k_size" : 10,
    "m" : 2,
    "n" : 3
}
```

#### 2.3.2 Extract Minimzers with k = 5, w = 10

```python
parse_params = {
    "mode" : 3,
    "k_size" : 5,
    "w" : 10,
}
```

