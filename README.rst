==========
KProcessor
==========

kProcessor is a C++ API with a handy python interface that enables easy handling of sequence kmer content. kProcessor stores kmers with their associated metadata in a virtual data structure called kDataframe. By default, kDataFrame stores the kmers with their counts in the input dataset. However, kDataFrame supports adding multiple other columns to store more information about these kmers in different data types. Users can easily merge multiple kDataFrames or apply different set functions (e.g. union, intersect and difference) on a group of kDataFrames. Colored kDataFrame is another core virtual data structure in kProcessor that allows the indexing of the kmers in a multi-sequence reference input. It is composed of a kDataFrame that replaces the kmer count with a key (aka color). This key connects the kmer to all sequences associated with this kmer.

**Types of kDataframe:**
kDataframe is a virtual data structure that has multiple implementations to fit different applications. Each implementation incorporates one of multiple data structures with different capabilities and limitations. Currently, kProcessor API has 4 versions of kDataframes. In the exact mode, all four versions should not be used to store kmers > 31 base pairs. 

*    **kDataFrameMQF** incorporates the mixed-counter quotient filter (`MQF <https://github.com/dib-lab/MQF>`_). MQF is an extremely fast and light-weight data structure for the kmers insertion and querying. The kmer hash values (not the kmers themselves) are sorted in MQF which enables easy implementation of set functions. The current implementation can only store kmers greater than 17 base pairs. 
*    **kDataFrameBMQF** incorporates BufferedMQF which is a version of MQF that trades some of the speed of insertions and queries for a significant memory reduction by scaling out of the memory to SSD disk
*    **kDataFramePHMAP** incorporates Parallel Hash Map. This unsorted data structure comes second in speed and memory efficiency after MQF but useful for small kmers. 
*    **kDataFrameMAP** incorporates MAP as its data structure. MAP is a balanced tree structure that allows the kmers to be sorted lexicographically. It can store small kmers as well. However, kDataFrameMAP is the most computationally intensive kDataframe regarding memory and time consumption. 

**Data input in kDataframe:**
kDataframes store the items as hash values which are unsigned integers. However they receive their inputs in different formats for convenience:

a.    **Unsigned integers:** This option allows the users to encode their input items in their own way and enables easy integration in different software packages. 
b.    **Kmers, sequences or sequence files:** Any kDataframe can parse ACGT sequences in different hashing modes (e.g. exact or approximate hashing), kmer orientations (canonical or non-canonical), and parsing modes (e.g. classical kmers, skipmers, or minimizers). Read more about these options at the kDataframe construction and sequence parsing functions.    


.. image:: https://travis-ci.org/dib-lab/kProcessor.svg?branch=master
    :target: https://travis-ci.org/dib-lab/kProcessor

.. .. image:: _static/images/kDataFrame.png
