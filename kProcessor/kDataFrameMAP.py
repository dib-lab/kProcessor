
from kProcessor.kDataFrame import kDataFrame

class kDataFrameMAP(kDataFrame):
    """
    The abstract base class defining a kDataFrameMAP.
    """

    def __init__(self, kSize):
        """
        :param kSize: Kmer Size
        :type kSize: integer
        :return: :class:`kProcessor.kDataFrameMAP`

        .. note:: Read more about the usage of kDataFrameMAP in the FAQ page.

        """
        pass

    def getTwin(self):
        """creates a new ``kDataFrameMAP`` using the same parameters as the current ``kDataFrameMAP``.

        :return: A shallow copy of the current ``kDataFrameMAP``.
        :rtype: kDataFrameMAP

        """
        pass

    def reserve(self, n):
        """Request a capacity change so that the kDataFrameMAP can approximately hold at least n kmers

        :param n: Minimum number of kmers
        :type n: integer
        """
        pass

    def insert(self, kmer, N = 1):
        """insert the kmer N time in the kDataFrameMAP, or increment the kmer count with N if it is already exists.

        :param kmer: The Kmer to increment its count
        :type kmer: string
        :param N: Kmer count (Optional, Default = 1)
        :type N: integer
        :return: Boolean value indicating whether the kmer is inserted or not
        :rtype: boolean
        """

        pass

    def setCount(self, kmer, N):
        """set the kmer's count to N time in the kDataFrameMAP

        :param kmer: The Kmer to set its count
        :type kmer: string
        :param N: Kmer count
        :type N: integer
        :return: Boolean value indicating whether the kmer is inserted or not
        :rtype: boolean
        """
        pass

    def count(self, kmer):
        """retrieve number of times the kmer was inserted in the kDataFrameMAP

        :param kmer: The kmer to retrieve its count
        :type kmer: string
        :return: The count of the kmer in the kDataFrameMAP
        :rtype: integer
        """
        pass

    def erase(self, kmer):
        """Removes  a kmer from the kDataFrameMAP

        :param kmer: The kmer to be erased
        :type kmer: string
        :return: Boolean value indicating whether the kmer is erased or not
        :rtype: boolean
        """
        pass

    def size(self):
        """ Number of kmers in the kDataFrameMAP

        :return: The number of kmers in the kDataFrameMAP
        :rtype: integer
        """
        pass

    def max_size(self):
        """ Maximum number of kmers that the kDataFrameMAP can hold.

        :return: The maximum number of kmers that the kDataFrameMAP can hold.
        :rtype: integer
        """
        pass

    def empty(self):
        """ Check whether the kDataFrameMAP is empty of kmers or not.

        :return: Boolean value indicating whether the kDataFrameMAP is empty, i.e. whether its size is 0
        """
        pass

    def load_factor(self):
        """ Retrieving the load factor of the kDataFrameMAP

        :return: The current load factor in the kDataFrameMAP.
        :rtype: integer
        """
        pass

    def max_load_factor(self):
        """ Retrieving the maximum load factor of the kDataFrameMAP

        :return: The maximum load factor in the kDataFrameMAP.
        :rtype: integer
        """
        pass

    def begin(self):
        """ Instantiate a :class:`kProcessor.kDataFrameIterator` object pointing to the first kmer position

        :return: An iterator at the begin of the kDataFrameMAP.
        :rtype: :class:`kProcessor.kDataFrameIterator`
        """

        pass

    def end(self):
        """ Instantiate a :class:`kProcessor.kDataFrameIterator` object pointing to the last kmer position

        :return: An iterator at the end of the kDataFrameMAP.
        :rtype: :class:`kProcessor.:class:`kProcessor.kDataFrameIterator``
        """
        pass

    def save(self):
        """
        Save the kDataFrameMAP on the disk in a form of binary file alongside other metadata files.
        Extension: ``.map``
        """
        pass

    def getKmerDecoder(self):
        """ Get the kmerDecoder instance object that's initialized in the kDataFrameMAP

        :return: The kmerDecoder instance used by kDataFrameMAP
        """
        pass

    @staticmethod
    def load(filePath):
        """ A static method to load a kDataFrameMAP file from disk.

        .. note:: Load the file without the extension [.map]

        :param filePath: The serialized kDataFrameMAP binary file without the extension
        :return: the loaded kDataFrameMAP from disk
        :rtype: :class:`kProcessor.kDataFrameMAP`

        Example:
            >>> import kProcessor as kp
            >>> # File path : "path/to/file.map"
            >>> KF = kp.kDataFrameMAP.load("path/to/file")
        """

        pass

    def kSize(self):
        """
        Get the kmer size of the kDataFrameMAP

        :return: kmer size
        :rtype: integer
        """
        pass

    def setkSize(self, k):
        """
        set the kmer size of the kDataFrameMAP
        """
        pass