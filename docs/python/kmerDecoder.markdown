# kmerDecoder

kmerDecoder is important to tell the kProcessor how to parse sequences.

It can be initialized to parse sequences files "FASTA/Q" or sequence strings.

## Modes

There are three modes for decoding the sequences substrings {kmers, skipmers, minimizers}

## Parameters

Mode parameters is a dictionary contains the parameters to initialize the kmerDecoder instance.

1. "kmers" mode parameters : `{"k_size" : k }` , where k is an integer.
2. "skipmers" mode parameters : `{"k_size" : k, "m" : m, "n" : n }` , where k,m,n are integers.
3. "minimizers" mode parameters : `{"k_size" : k, "w" : w }` , where k,w are integers.

[read more about skipmers](https://www.biorxiv.org/content/biorxiv/early/2017/09/19/179960.full.pdf)

## Initialization

### Initialize to parse kmers

```python

# FASTA/Q File to be parsed
filename = "sample.fa"

kmer_size = 21

# The chunk size is to determine how many sequences to read and process at once
chunk_size = 1000

mode = "kmers"
params = {"k_size": kmer_size}

KD = kp.initialize_kmerDecoder(filename, chunk_size, mode, params)

```

### Initialize to parse single sequence string

```python

kmer_size = 21
mode = "kmers"
params = {"k_size": kmer_size}
KD = kp.initialize_kmerDecoder(mode, params)

```
