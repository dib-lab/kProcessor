==========
Python API
==========

The main programmatic way to interact with the ``kProcessor`` is via its Python API.
Please also see `recipes of using the python API <recipes.html>`_.

.. contents::
   :depth: 3

``kDataFrame``: The abstract base class for defining the kDataFrame
===================================================================

.. automodule:: kProcessor.kDataFrame
   :members:
   :private-members:


``kDataFrameIterator``: The abstract base class for defining a kDataFrame iterator
==================================================================================

.. automodule:: kProcessor.kDataFrameIterator
   :members:
   :private-members:

``kDataFrameMQF``: subclass derived from ``kDataFrame``
=======================================================

.. automodule:: kProcessor.kDataFrameMQF
   :members:

``kDataFrameMAP``: subclass derived from ``kDataFrame``
=======================================================

.. automodule:: kProcessor.kDataFrameMAP
   :members:


``kDataFramePHMAP``: subclass derived from ``kDataFrame``
=========================================================

.. automodule:: kProcessor.kDataFramePHMAP
   :members:


``colored_kDataFrame``: colored kDataFrame that holds the source sequence of each Kmer
======================================================================================

.. automodule:: kProcessor.colored_kDataFrame
   :members:


``Set Functions``: Function like intersection & union
=====================================================

.. automodule:: kProcessor.set_function
   :members:


``Sequence parsing``: Reading sequence kmers into kDataFrames
=============================================================

.. automodule:: kProcessor.seqs_to_kmers
   :members:
