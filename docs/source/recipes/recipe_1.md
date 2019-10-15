# Recipe 1: Kmers parsing and counting

## Data

Sample light-weight data for running the examples.

<a href="../../../assets/files/data.zip" download> Click here to download </a>

## Description

1. Create an empty kDataFrame with kmerSize = 21
2. Load a fasta file into a kDataFrame
3. Save the kDataFrame on disk

## Implementation

### Importing

```python
import kProcessor as kp
```

### Create an empty kDataFrame

```python
kf1 = kp.kDataFrameMQF(21)
```

### Create KmerDecoder Object
> MUST be initialized with the same kDataFrame kSize

```python
KD = kp.initialize_kmerDecoder("data/test.fastq" , 1000,"kmers", {"k_size" : 21})
```

### Parse the fastq file into the kf1 kDataFrame

```python
kp.parseSequences(KD, kf1)
```

### Iterating over first 10 kmers

```python
it = kf1.begin()

for i in range(10):
    print(it.getKmer())
    it.next() # Extremely important to move the iterator to the next kmer position.
```

### Save the kDataFrame on disk with a name "kf1"

```python
# This will save the file with the extension ".mqf"
kf1.save("kf1")
```

## Complete Script

```python
import kProcessor as kp

# Creating an empty kDataFrameMQF with kmer size 21
kf1 = kp.kDataFrameMQF(21)

KD = kp.initialize_kmerDecoder("data/test.fastq" , 1000,"kmers", {"k_size" : 21})

KD = kp.initialize_kmerDecoder("data/test.fastq" , 1000,"kmers", {"k_size" : 21})
kp.parseSequences(KD, kf1)
kf1.save("kf1")

it = kf1.begin()

for i in range(10):
    print(it.getKmer())
    it.next() # Extremely important to move the iterator to the next kmer position.

```

**Output**

```text
CCCAACAGAATTAAAAAGTCA
AAATTAAATAACTTTAGCGCA
CCAAATTACAACAAAATTTGG
TTAATCATTTGGTATAATTGC
ACCTCGTATAACTTCGTATAA
AACAATTCAACAGAGAAGGAC
AGGCTAATCGAACAAAACATC
AGGAAAAACTCCAGCCAGTAA
TACGGGTCGCAGTGACCAGGC
CCAGGTAGTACAGCAATCGTA
```