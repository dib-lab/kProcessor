# kDataFrame

`kp.kDataFrame` should be replaced with any instance of the type `kDataFrame` either it's type is `kp.kDataFrameMQF` or `kp.kDataFrameMAP`.

Example:

```python

import kProcessor as kp

kf = kDataFrameMQF(31)

# Get kDataFrame size

kf_size = kf.size()

```

## Notes

- The maximum **kmer size** can be handled is `31`
- `kDataFrameMQF` is faster than `kDataFrameMAP` and has lower memory consumption.
- `kDataFrameMAP` is used when processing with **kmer size** < `15`
- kDataFrame set functions could be applied only on `kDataFrameMQF`
- Each of the two subclasses provides the same functions.

---

### Instantiation

Create an empty kDataFrame with predefined kmerSize

```python

empty_kFrame = kp.kDataFrame(kSize)

```

---

## kDataFrame basic functions

### Get the size of the kDataFrame

returns the number of kmers in the kDataFrame

```python

kDataFrame.size()

```

### Get the maximum size

returns the maximum number of kmers that the kDataframe can hold.

```python

kDataFrame.max_size()

```

### Saving the kDataFrame on disk

In case of `kDataFrameMQF` will write on disk a file with the extension `.mqf`
In case of `kDataFrameMAP` will write on disk a file with the extension `.phmap`

```python

kDataFrame.save(fileName)

```

### Loading the kDataFrame from disk

- Load the kDataFrame From disk and return kDataFrame object with the same loaded kDataFrame type.
- The function `load(filename)` is static function,  don't replace `kDataFrame` class with `kDataFrameMQF` or `kDataFrameMAP`.
- Don't write the file extension `.mqf` or `.pham`, this is autodetected by the software.

```python

loaded_kf = kp.kDataFrame.load(filename)

```

---

## kDataFrame kmer-specific functions


kDataFrame.count(kmer): returns the kmer count (number of times that kmer was inserted)
kDataFrame.erase(kmer): remove a kmer from the kDataFrame

### Insert a kmer

Insert a kmer in the `kDataFrame`

```python

kDataFrame.insert(kmer)

```

### Insert a kmer and set its count

Insert a kmer alongside its count in the kDataFrame, this will **increment** the old count with the new one if the kmer is already there.

```python

kDataFrame.insert(kmer, count)

```

### set a kmer count

Insert a kmer alongside its count in the kDataFrame, this will **reset** the old count with the new one if the kmer is already there.

```python

kDataFrame.setCount(kmer, count)

```

### get the kmer count

returns the kmer count (number of times that kmer was inserted)

```python

kDataFrame.count(kmer)

```

### remove a kmer from the kDataFrame

remove a kmer from the kDataFrame

```python

kDataFrame.erase(kmer)

```
