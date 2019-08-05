# kmer counting

## kmer counting from sequence string

```python

kp.parseSequencesFromString(kmerDecoder, sequence_string, kDataFrame)

```

### Example

```python

# 1- Instantiate kmerDecoder object
KD = kp.initialize_kmerDecoder("kmers", {"k_size":21})

# 2 - Instantiate kDataFrame
KF = kp.kDataFrameMQF(21)

# 3- Parse the sequene and insert the kmers in the kDataFrame
kp.parseSequencesFromString(KD, seq, KF)


```

## kmer counting from FASTA/Q file

```python

kp.parseSequences(kmerDecoder, kDataFrame)

```

### Example

```python

# 1- Instantiate kmerDecoder object
KD = kp.initialize_kmerDecoder("samples.fa", 1000, "kmers", {"k_size":21})

# 2 - Instantiate kDataFrame
KF = kp.kDataFrameMQF(21)

# 3- Parse the sequenes and insert the kmers in the kDataFrame
kp.parseSequencesFromString(KD, KF)

```
