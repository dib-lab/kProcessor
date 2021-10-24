# FAQ

<hr>

!!! note
    Hashing and Reading, modes alongside with the kmers parameters are tightly associated with the kDataFrame object.

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

## **2. Sequence Parsing**

In order to extract the kmers from a sequences file or string in kProcessor, you will need to pass a Python dictionary with the parsing parameters.

### 2.1 Parameters

Keys of the dictionary are *strings* describing the parameter type, values are integers for the parameters values

### 2.2 Modes

There are four modes for decoding the sequences substrings `KMERS`, `SKIPMERS`, `MINIMIZERS`, and `PROTEIN`. The following table describes the reading mode/hashing mode compatibility. True values represent the modes that can be used together.

| **Reading Mode** | **Hashing Mode**           | **Compatible?** |
|:----------------:|:--------------------------:|:---------------:|
| KMERS            | mumur_hasher               | true            |
| SKIPMERS         | mumur_hasher               | true            |
| MINIMIZERS       | mumur_hasher               | true            |
| PROTEIN          | mumur_hasher               | false           |
| KMERS            | integer_hasher             | true            |
| SKIPMERS         | integer_hasher             | true            |
| MINIMIZERS       | integer_hasher             | true            |
| PROTEIN          | integer_hasher             | false           |
| KMERS            | TwoBits_hasher             | true            |
| SKIPMERS         | TwoBits_hasher             | true            |
| MINIMIZERS       | TwoBits_hasher             | true            |
| PROTEIN          | TwoBits_hasher             | false           |
| KMERS            | nonCanonicalInteger_Hasher | true            |
| SKIPMERS         | nonCanonicalInteger_Hasher | true            |
| MINIMIZERS       | nonCanonicalInteger_Hasher | true            |
| PROTEIN          | nonCanonicalInteger_Hasher | false           |
| KMERS            | protein_hasher             | false           |
| SKIPMERS         | protein_hasher             | false           |
| MINIMIZERS       | protein_hasher             | false           |
| PROTEIN          | protein_hasher             | true            |
| KMERS            | proteinDayhoff_hasher      | false           |
| SKIPMERS         | proteinDayhoff_hasher      | false           |
| MINIMIZERS       | proteinDayhoff_hasher      | false           |
| PROTEIN          | proteinDayhoff_hasher      | true            |

#### 2.2.1 KMERS
    
- Description: Extracts the sequences substrings in the default popular mode "kmers".
- Reading Mode: `KMERS`
- Parameters: 
    - k_size: total number of nucleotide bases

#### 2.2.2 SKIPMERS

- Description: A cyclic pattern of picked or skipped positions. [read more about skipmers](https://www.biorxiv.org/content/10.1101/179960v2)
- Reading Mode: `SKIPMERS`
- Parameters:
    - "k_size": total number of bases 
    - "m": used bases per cycle
    - "n": cycle length

#### 2.2.3 MINIMIZERS

- Description: short substrings that represents the sequence. [read more about minimizers](https://homolog.us/blogs/bioinfo/2017/10/25/intro-minimizer/)
- Reading Mode: `MINIMIZERS`
- Parameters:
    - "k_size": total number of bases 
    - "w": window size

#### 2.2.4 PROTEIN

- Description: Parsing protein sequences to extract the sequence substrings.
- Reading Mode: `PROTEIN`
- Parameters: 
    - k_size: total number of amino acid bases

### 2.3 Default Hashing and Reading Modes

- kDataFrameMQF/kDataFrameBMQF: Default reading mode is `KMERS` and the default hashing mode is `integer_hasher` which is reversible.
- kDataFramePHMAP: Default reading mode is `KMERS` and the default hashing mode is `TwoBits_hasher` which is reversible.
- kDataFrameMAP: Default reading mode is `KMERS` and the default hashing mode is `TwoBits_hasher` which is reversible.

!!! note "technical information"
    The `twoBits_hasher` is used to avoid the double hashing of the kmers since the underlying data structures has an internal hashing. `twoBits_hasher` is not an actual hasher, it just converts a kmer substring to it's corresponding two-bits representation.


### 2.4 Examples

#### 2.4.1 Extract kmers with kmer size 31

```python
parse_params = {
    "kSize" : 31
}

KF_KMERS = kDataFramePHMAP(KMERS, twoBits_hasher, parse_params)
```

#### 2.4.2 Extract skipmers with k = 10, m = 2, n = 3

```python
parse_params = {
    "k" : 10,
    "m" : 2,
    "n" : 3
}

KF_SKIPMERS = kp.kDataFramePHMAP(SKIPMERS, integer_hasher, parse_params)
```

#### 2.4.3 Extract Minimzers with k = 5, w = 10

```python
parse_params = {
    "k" : 5,
    "w" : 10,
}

KF_MINIMIZERS = kp.kDataFramePHMAP(MINIMIZERS, twoBits_hasher, parse_params)
```

#### 2.4.4 Extract Protein with kSize = 5

```python
parse_params = {
    "kSize" : 5,
}

KF_PROTEIN = kp.kDataFramePHMAP(PROTEIN, protein_hasher, parse_params)
```
