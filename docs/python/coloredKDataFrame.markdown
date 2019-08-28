# coloredKDataFrame

## Introduction

- `coloredKDataFrame` class inherits all the functions from `kDataFrame` plus other new functions.
- The `coloredKDataFrame` class holds the kmers colors instead of their count. so, `coloredKDataFrame.count(kmer)` will return the color.
- A kmer color is an integer represents the targets which contains that kmer.

Example:

`color : 1`: represents the transcripts `transcript_A` , `transcript_B` and `transcript_C`
`color : 2`: represents the transcripts `transcript_A` , `transcript_B`

`kmer: ACTGATCGATCGTACGAC` has the color `2`, that means it's found in both `transcript_A` and `transcript_B`
`kmer: ATAAGCATTTACAGCAAT` has the color `1`, that means it's found in both `transcript_A` , `transcript_B` and `transcript_C`

---


## Functions

### Get kmer color

returns the color of the kmer

```python

kmer_color = colored_kDataFrame.getKmerColor(kmer)

```

### Get all samples (targets) ID(s) contains a given color

return a list of samples IDs associated with this color

```python

samples_list = colored_kDataFrame.getSamplesIDForColor(color)

```

### Get all colors of a given kmer

return a list of colors associated with this kmer

```python

colors_list = colored_kDataFrame.getSamplesIDForKmer(kmer)

```

### Get the names map

returns a dictionary of names map

```python

colored_kDataFrame.names_map()

```

### Get the inversed names map

returns an inversed dictionary of names map

```python

colored_kDataFrame.inversed_names_map()

```
