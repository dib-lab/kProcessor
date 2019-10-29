# FAQ

<hr>

## **Hashing**

### What are kDataFrame hashing modes?

> **Hashing modes are only applied in kDataFrameMQF and kDataFramePHMAP**

There are multiple hasing modes for the `kDataFrame` to use.
- 0: When setting the mode parameter to `0` that will define the hashing to be irreversible. By other words, once the kmer is hashed and inserted in the `kDataFrame` can't be reversed back to it's string representation.
- 1: Hashing mode `1` will define the hashing to be reversible, so that you can iterate over the `kDataFrame` kmers in their string representation. 

### Why there are not user-selected hashing modes in `kDataFrameMAP`?

The only purpose of using the `kDataFrameMAP` is to store the kmers lexicographically. Consequently, no hashing is required, just store the kmers in their Two-bits representation.

### What's the best hashing mode to use?

There's no best hashing mode, that depends totally on the application.

The irreversible mode is way faster than the reversible mode and requires less space. However, it's not as sensitive as the reversible mode in hashing the kmers because it stores the approximate hash value of the kmers with a user-controlled false-positive rate.
As the false-positive rate increases, the probability of hashing multiple kmers with the same hash value will increase.
Moreover, querying the kDataFrame with a similar but not stored kmer can show that the kmer already exist, but actually it's not. 

On the other hand, the reversible hashing mode is slower than the irreversible mode but allows the user to iterate over the kmers in their string representation.

---

## **Section two**

### Here goes question one
> Answer of question one

---

## **Section three**

### Here goes question one
> Answer of question one

