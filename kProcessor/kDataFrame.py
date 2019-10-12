"""
.. module:: kProcessor
"""

class kDataFrame:
    """
    The abstract base class defining a kDataFrame.
    """

    def __init__(kSize):
        """kDataFrame constructor.

        Instantiate a kDataFrame object with predefined kmer size.

        :param kSize: Kmer Size
        :type kSize: integer
        :return: Void

        """
        pass

    def getTwin(self):
        """creates a new ``kDataframe`` using the same parameters as the current ``kDataFrame``.

        :return: A shallow copy of the current ``kDataFrame``.
        :rtype: kDataFrame

        """
        pass

    def reserve(self, n):
        """Request a capacity change so that the kDataFrame can approximately hold at least n kmers

        :param n: Minimum number of kmers
        :type n: integer
        """
        pass

    def insert(self, kmer, N = 1):
        """insert the kmer N time in the kDataFrame, or increment the kmer count with N if it is already exists.

        :param kmer: The Kmer to increment its count
        :type kmer: string
        :param N: Kmer count (Optional, Default = 1)
        :type N: integer
        :return: Boolean value indicating whether the kmer is inserted or not
        :rtype: bool
        """

        pass

    def setCount(self, kmer, N):
        """set the kmer's count to N time in the kDataFrame

        :param kmer: The Kmer to set its count
        :type kmer: string
        :param N: Kmer count
        :type N: integer
        :return: Boolean value indicating whether the kmer is inserted or not
        :rtype: bool
        """
        pass

    def count(self, kmer):
        """retrieve number of times the kmer was inserted in the kdataFrame

        :param kmer: The kmer to retrieve its count
        :type kmer: string
        :return: The count of the kmer in the kDataFrame
        :rtype: integer
        """
        pass

    def erase(self, kmer):
        """Removes  a kmer from the kDataFrame

        :param kmer: The kmer to be erased
        :type kmer: string
        :return: Boolean value indicating whether the kmer is erased or not
        :rtype: bool
        """
        pass

    def size(self):
        """ Number of kmers in the kDataFrame

        :return: The number of kmers in the kDataframe
        :rtype: integer
        """
        pass

    def max_size(self):
        """ Maximum number of kmers that the kDataframe can hold.

        :return: The maximum number of kmers that the kDataframe can hold.
        :rtype: integer
        """
        pass

    def empty(self):
        """ Check whether the kDataFrame is empty of kmers or not.

        :return: Boolean value indicating whether the kDataFrame is empty, i.e. whether its size is 0
        """
        pass

    def load_factor(self):
        """ Retrieving the load factor of the kDataFrame

        :return: The current load factor in the kDataFrame.
        :rtype: integer
        """
        pass

    def max_load_factor(self):
        """ Retrieving the maximum load factor of the kDataFrame

        :return: The maximum load factor in the kDataFrame.
        :rtype: integer
        """
        pass

    def begin(self):
        """ Instantiate a kDataFrameIterator object pointing to the first kmer position
        :return: An iterator at the begin of the kDataFrame.
        :rtype: :class:`kProcessor.kDataFrameIterator`
        """

        pass

    def end(self):
        """ Instantiate a kDataFrameIterator object pointing to the last kmer position

        :return: An iterator at the end of the kDataFrame.
        :rtype: :class:`kProcessor.kDataFrameIterator`

        """
        pass

    def save(self):
        """
        Save the kDataFrame on the disk in a form of binary file alongside other metadata files.
        """
        pass

    def getKmerDecoder(self):
        """ Get the kmerDecoder instance object that's initialized in the kDataFrame

        :return: The kmerDecoder instance used by kDataframe
        :rtype: kmerDecoder
        """
        pass

    @staticmethod
    def load(filePath):
        """ A static method to load a kDataFrame file from disk.

        .. note:: Load the file without the extension [.mqf, .map, .phmap]

        :param filePath: The serialized kDataFrame binary file without the extension
        :return: the loaded kDataFrame from disk
        :rtype: :class:`kProcessor.kDataFrame`

        Example:
            >>> import kProcessor as kp
            >>> # File path : "path/to/file.mqf"
            >>> KF = kp.kDataFrame.load("path/to/file")
        """

        pass

    def kSize(self):
        """
        Get the kmer size of the kDataFrame

        :return: kmer size
        :rtype: integer
        """
        pass

    def setkSize(self, k):
        """
        set the kmer size of the kDataFrame
        """
        pass