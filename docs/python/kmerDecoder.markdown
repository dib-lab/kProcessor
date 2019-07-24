# kmerDecoder

kmerDecoder is important to tell the kProcessor how to parse sequences.

It can be initialized to parse sequences files "FASTA/Q" or sequence strings.

## Modes

There are three modes for decoding the kmers {kmers, skipmers, minimizers}

## Parameters

mode parameters is a dictionary contains the parameters to initialize the kmerDecoder instance.

"kmers" mode parameters : `{"k_size" : k }` , where k is an integer.
"skipmers" mode parameters : `{"k_size" : k, "m" : m, "n" : n }` , where k,m,n are integers.
"minimizers" mode parameters : `{"k_size" : k, "w" : w }` , where k,w are integers.

[read more about skipmers](https://www.biorxiv.org/content/biorxiv/early/2017/09/19/179960.full.pdf)

## Initialization example

```python

filename = "sample.fa"
kmer_size = 21
mode = "kmer"

kp.initialize_kmerDecoder()

```
