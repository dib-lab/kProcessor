
from kProcessor.kDataFrame import kDataFrame

class kDataFrameMQF(kDataFrame):
    """The abstract base class defining a kDataFrameMQF.
    """

    def __init__(self, kSize, mode=1):
        """Instantiate a kDataFrameMQF object with predefined kmer size.

        :param kSize: Kmer Size
        :type kSize: integer
        :param mode: Hashing mode for the kDataFrameMQF, default  = 1
        :type mode: integer
        :return: :class:`kProcessor.kDataFrameMQF`

        Instantiation Example:
            >>> import kProcessor as kp
            >>> KF_MQF = kp.kDataFrameMQF(31) # kSize = 31, Hashing Mode = 1
            >>> KF2_MQF = kp.kDataFrameMQF(31, 0) # kSize = 31, Hashing Mode = 0

        .. note:: Read more about hashing modes in the FAQ page.

        """
        pass


    def getTwin(self):
        """creates a new ``kDataFrameMQF`` using the same parameters as the current ``kDataFrameMQF``.

        :return: A shallow copy of the current ``kDataFrameMQF``.
        :rtype: kDataFrameMQF

        """
        pass

    def reserve(self, n):
        """Request a capacity change so that the kDataFrameMQF can approximately hold at least n kmers

        :param n: Minimum number of kmers
        :type n: integer
        """
        pass

    def insert(self, kmer, N = 1):
        """insert the kmer N time in the kDataFrameMQF, or increment the kmer count with N if it is already exists.

        :param kmer: The Kmer to increment its count
        :type kmer: string
        :param N: Kmer count (Optional, Default = 1)
        :type N: integer
        :return: Boolean value indicating whether the kmer is inserted or not
        :rtype: boolean
        """

        pass

    def setCount(self, kmer, N):
        """set the kmer's count to N time in the kDataFrameMQF

        :param kmer: The Kmer to set its count
        :type kmer: string
        :param N: Kmer count
        :type N: integer
        :return: Boolean value indicating whether the kmer is inserted or not
        :rtype: boolean
        """
        pass

    def count(self, kmer):
        """retrieve number of times the kmer was inserted in the kDataFrameMQF

        :param kmer: The kmer to retrieve its count
        :type kmer: string
        :return: The count of the kmer in the kDataFrameMQF
        :rtype: integer
        """
        pass

    def erase(self, kmer):
        """Removes  a kmer from the kDataFrameMQF

        :param kmer: The kmer to be erased
        :type kmer: string
        :return: Boolean value indicating whether the kmer is erased or not
        :rtype: boolean
        """
        pass

    def size(self):
        """ Number of kmers in the kDataFrameMQF

        :return: The number of kmers in the kDataFrameMQF
        :rtype: integer
        """
        pass

    def max_size(self):
        """ Maximum number of kmers that the kDataFrameMQF can hold.

        :return: The maximum number of kmers that the kDataFrameMQF can hold.
        :rtype: integer
        """
        pass

    def empty(self):
        """ Check whether the kDataFrameMQF is empty of kmers or not.

        :return: Boolean value indicating whether the kDataFrameMQF is empty, i.e. whether its size is 0
        """
        pass

    def load_factor(self):
        """ Retrieving the load factor of the kDataFrameMQF

        :return: The current load factor in the kDataFrameMQF.
        :rtype: integer
        """
        pass

    def max_load_factor(self):
        """ Retrieving the maximum load factor of the kDataFrameMQF

        :return: The maximum load factor in the kDataFrameMQF.
        :rtype: integer
        """
        pass

    def begin(self):
        """ Instantiate a :class:`kProcessor.kDataFrameIterator` object pointing to the first kmer position

        :return: An iterator at the begin of the kDataFrameMQF.
        :rtype: :class:`kProcessor.kDataFrameIterator`
        """

        pass

    def end(self):
        """ Instantiate a :class:`kProcessor.kDataFrameIterator` object pointing to the last kmer position

        :return: An iterator at the end of the kDataFrameMQF.
        :rtype: :class:`kProcessor.:class:`kProcessor.kDataFrameIterator``
        """
        pass

    def save(self):
        """
        Save the kDataFrameMQF on the disk in a form of binary file alongside other metadata files.
        Extension: ``.mqf``
        """
        pass

    def getKmerDecoder(self):
        """ Get the kmerDecoder instance object that's initialized in the kDataFrameMQF

        :return: The kmerDecoder instance used by kDataFrameMQF
        """
        pass

    @staticmethod
    def load(filePath):
        """ A static method to load a kDataFrameMQF file from disk.

        .. note:: Load the file without the extension [.mqf]

        :param filePath: The serialized kDataFrameMQF binary file without the extension
        :return: the loaded kDataFrameMQF from disk
        :rtype: :class:`kProcessor.kDataFrameMQF`

        Example:
            >>> import kProcessor as kp
            >>> # File path : "path/to/file.mqf"
            >>> KF = kp.kDataFrameMQF.load("path/to/file")
        """

        pass

    def kSize(self):
        """
        Get the kmer size of the kDataFrameMQF

        :return: kmer size
        :rtype: integer
        """
        pass

    def setkSize(self, k):
        """
        set the kmer size of the kDataFrameMQF
        """
        pass