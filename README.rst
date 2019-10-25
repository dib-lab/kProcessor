==========
KProcessor
==========

kProcessor is a C++ API with a handy python interface that enables easy handling of sequence kmer content. kProcessor stores kmers with their associated metadata in a virtual data structure called kDataframe. By default, kDataFrame stores the kmers with their counts in the input dataset. However, kDataFrame supports adding multiple other columns to store more information about these kmers in different data types. Users can easily merge multiple kDataFrames or apply different set functions (e.g. union, intersect and difference) on a group of kDataFrames. Colored kDataFrame is another core virtual data structure in kProcessor that allows the indexing of the kmers in a multi-sequence reference input. It is composed of a kDataFrame that replaces the kmer count with a key (aka color). This key connects the kmer to all sequences associated with this kmer.

kDataframe has multiple implementations to fit different applications:
======================================================================
kDataframe is a virtual data structure that incorporates one of multiple data structures, each one of these structures has different capabilities and limitations. Currently, the kProcessor API has 4 versions of the kDataframe. The kDataFrameMQF uses the mixed-counter quotient filter (MQF). MQF is an extremely fast and light-weight data structure for the kmers insertion and querying. The kmer hash values (not the kmers themselves) are sorted in MQF which enables easy implementation of set functions. The current implementation can store kmers greater than 17 base pairs. The kDataFramePHMAP derived class uses the Parallel Hash Map. This unsorted data structure comes second in speed and memory efficiency after MQF but useful for small kmers. Finally, The kDataFrameMAP is another kDataframe that uses MAP as its data structure. MAP is a balanced tree structure that allows the kmers to be sorted lexicographically. It can store small or even mixed size kmers. However, kDataFrameMAP is the most computationally intensive kDataframe regarding memory and time consumption. In the exact mode, all four versions should not be used to store kmers > 31 base pairs. 


Data represenation in kDataframe: 
=================================
kDataframes store the items as hash values which are unsigned integers. Users can insert integers directly into a kDataframe however the typical inputs are kmers, sequences or sequence files. Any kDataframe has an immutable hashing-mode which allows the storage of kmers in either an exact reversible or an approximate irreversible representations. The approximate mode is more memory efficient but comes with the risk of a small error rate. They kmers are typically stored in the canonical format where a kmer and its reverse complement have the same hash value, however the user can change this setting the kDataframe in the non-canonical mode. kProcessor has different parsing-modes of sequences. The classical parsing-mode is designed to convert each given sequence into a set of overlapping kmers where each two sequential kmers share k-1 base pairs. The skipper mode generate similar kmers but the parser skip n bases every m bases (n and m are user defined). The minimizer mode generate a representative kmer for each window w of sequences (w is user defined).


.. image:: https://travis-ci.org/dib-lab/kProcessor.svg?branch=master
    :target: https://travis-ci.org/dib-lab/kProcessor

.. image:: _static/images/kDataFrame.png
