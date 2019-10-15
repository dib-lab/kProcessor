# Recipe 2: Loading a kDataFrame file from disk.

## Description

In [Recipe 1](./recipe_1.md) we saved a kDataFrame to the disk with the file name `kf1.mqf`

1. Load the kf1 kDataFrame in an new kDataFrame.
2. Verify the loading by printing the kmerSize and total size.

## Implementation

### Importing

```python
import kProcessor as kp
```

### Loading the kDataFrame

```python
kf = kp.kDataFrame.load("kf1") # Note: We didn't write the extension.
```

### Print the kmer size

```python
kSize = kf.kSize()

print(f"kSize: {kSize}")
```

### Dump the first 10 kmers

```python
it = kf.begin()

for i in range(10):
    print(it.getKmer())
    it.next() # Extremely important to move the iterator to the next kmer position.
```

## Complete Script

```python
import kProcessor as kp

kf = kp.kDataFrame.load("kf1") # Note: We didn't write the extension.

kSize = kf.ksize()

print(f"kSize: {kSize}")

it = kf.begin()

for i in range(10):
    print(it.getKmer())
    it.next() # Extremely important to move the iterator to the next kmer position.

```

**Output**

```text
kSize: 21
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