# kDataFrame Iterators

kDataFrameIterator class is resposible for iterating over any kDataFrame instance to get pairs of kmer and its count in case of the `kDataFrameMAP` && `kDataFrameMQF`, the count will be replaced with color in the case of `coloredKDataFrame`.

## Iterator functions

### First position

Return a pointer to the first position of a kDataFrame object

```python

kDataFrameIterator.begin()

```

### Last position

Return a pointer to the last position of a kDataFrame object

```python

kDataFrameIterator.end()

```

### Move to next position

Increment the iterator position by one

```python

kDataFrameIterator.next()

```

### Get the kmer of the current iterator position

Return the kmer at the current iterator position

```python

kDataFrameIterator.getKmer()

```

### Get the count (or color in case of `coloredKDataFrame` ) of the current iterator position

Return the kmer at the current iterator position

```python

kDataFrameIterator.getCount()

```

### Set the kmer count (or color in case of `coloredKDataFrame` ) of the current iterator position

Set the kmer count at the current iterator position

```python

kDataFrameIterator.setCount(new_count)

```
