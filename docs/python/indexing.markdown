# Indexing

## Important terms

### kmer color

A color is a term used to represent all the targets *(indexed sequences)* that associated with a kmer.

#### Example

- kmer `AATTTT` found in `seq1`
- kmer `GGGGGG` found in `seq2`
- kmer `GCTGGG` found in `seq3`
- kmer `AACCGT` found in `seq1` & `seq2`
- kmer `GGTTAA` found in `seq2` & `seq3`

**After indexing these will be the colors**

- color `1` is associated with `seq1`
- color `2` is associated with `seq2`
- color `3` is associated with `seq3`
- color `4` is associated with `seq1` & `seq2`
- color `5` is associated with `seq2` & `seq3`

**That means:**
- color of the kmer `AATTTT` is **`1`** 
- color of the kmer `GGGGGG` is **`2`**
- color of the kmer `GGGGGG` is **`3`**
- color of the kmer `AACCGT` is **`4`**
- color of the kmer `GGTTAA` is **`5`**

---

### Names File

Names file is **tab delimited file (TSV)** passed alongside with the sequences file to demonstrate how the index will be held.

First column represents the exact FASTA header for each sequence in the without the `>` character.
Second column represents the corresponding group to that sequence.

#### Example 1

In the following names map file example, each transcript will be indexed a unique group.

```tsv

transcript1	transcript1
transcript2	transcript2
transcript3	transcript3
transcript4	transcript4

```

#### Example 2

In the following names map file example, `transcript1` and `transcript2` will be indexed as a single group called `group1` so all the kmers found with the `transcript1` and `transcript2` will be associated with `group1`.

```tsv

transcript1	group1
transcript2	group1
transcript3	group2
transcript4	group2

```

---

## Sequences Indexing

`kp.index()` takes a kmerDecoder object, names file path and a kDataFrame as input parameters. Returns a colored kDataFrame.

```python

ckf = kp.index(kmerDecoder, names_file_path, kDataFrame)

```

### Indexing Example

```python

fasta_file = "sample.fa"

names_file = "sample.fa.names"

kmer_size = 21

mode = "kmers"

params = {"k_size": 31}

chunk_size = 10000

# 1- Instantiate kmerDecoder object
KD = kp.initialize_kmerDecoder(fasta_file, chunk_size, "kmers", {"k_size" : kmer_size})

# 2 - Instantiate kDataFrame
KF = kp.kDataFrameMQF(kmer_size)

# 3- Perform indexing, store the kmers with their colors in the KF object and the colors information will be returned as a colored_kDataFrame (cfk)
cfk = kp.index(KD, names_file, KF)

# 4- Save the index to the disk
cfk.save("idx_sample")

```
