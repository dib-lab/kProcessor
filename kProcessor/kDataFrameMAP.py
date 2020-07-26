
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

    # Will publish it later after cloning kmerDecoder settings.
    # def getTwin(self):
    #     """creates a new ``kDataFrameMAP`` using the same parameters as the current ``kDataFrameMAP``.
    #
    #     :return: A shallow copy of the current ``kDataFrameMAP``.
    #     :rtype: kDataFrameMAP
    #
    #     """
    #     # pass

    def reserve(self, n):
        """Request a capacity change so that the kDataFrameMAP can approximately hold at least n kmers

        :param n: Minimum number of kmers
        :type n: integer

         .. note:: Read more about the usage of `reserve(n)` in the FAQ page.

        """
        pass

    def insert(self, kmer, N = 1):
        """Insert the kmer N time in the kDataFrameMAP, or increment the kmer count with N if it is already exists.

        :param kmer: The Kmer to increment its count
        :type kmer: string
        :param N: Kmer count (Optional, Default = 1)
        :type N: integer
        :return: Boolean value indicating whether the kmer is inserted or not
        :rtype: bool
        """

        pass

    def setCount(self, kmer, N):
        """Set the kmer's count to N time in the kDataFrameMAP

        :param kmer: The Kmer to set its count
        :type kmer: string
        :param N: Kmer count
        :type N: integer
        :return: Boolean value indicating whether the kmer is inserted or not
        :rtype: bool
        """
        pass

    def getCount(self, kmer):
        """Retrieve number of times the kmer was inserted in the kDataFrameMAP

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
        :rtype: bool
        """
        pass

    def size(self):
        """Number of kmers in the kDataFrameMAP

        :return: The number of kmers in the kDataFrameMAP
        :rtype: integer
        """
        pass

    def max_size(self):
        """Maximum number of kmers that the kDataFrameMAP can hold.

        :return: The maximum number of kmers that the kDataFrameMAP can hold.
        :rtype: integer
        """

        pass

    def empty(self):
        """Check whether the kDataFrameMAP is empty of kmers or not.

        :return: Boolean value indicating whether the kDataFrameMAP is empty, i.e. whether its size is 0
        :rtype: boolean

        """
        pass

    def load_factor(self):
        """Retrieving the current load factor of the kDataFrameMAP in percentage to indicate how full is it.

        :return: The current load factor in the kDataFrameMAP.
        :rtype: integer
        """
        pass

    def max_load_factor(self):
        """ Retrieving the maximum load factor of the kDataFrameMAP in percentage.

        :return: The maximum load factor in the kDataFrameMAP.
        :rtype: integer
        """
        pass

    def begin(self):
        """ Instantiate a kDataFrameIterator object pointing to the first kmer position
        :return: An iterator at the begin of the kDataFrameMAP.
        :rtype: :class:`kProcessor.kDataFrameIterator`
        """

        pass

    def end(self):
        """ Instantiate a kDataFrameIterator object pointing to the last kmer position

        :return: An iterator at the end of the kDataFrameMAP.
        :rtype: :class:`kProcessor.kDataFrameIterator`

        """
        pass

    def save(self):
        """
        Serialize the kDataFrameMAP on the disk in a form of binary file alongside other metadata files.
        """

        pass


    @staticmethod
    def load(filePath):
        """ A static method to load a kDataFrameMAP file from disk.

        .. note:: Load the file without the extension [.mqf, .map, .phmap]

        :param filePath: The serialized kDataFrameMAP binary file without the extension
        :return: the loaded kDataFrameMAP from disk
        :rtype: :class:`kProcessor.kDataFrameMAP`

        Example:
            >>> import kProcessor as kp
            >>> # File path : "path/to/file.mqf"
            >>> KF = kp.kDataFrameMAP.load("path/to/file")
        """

        pass

    def kSize(self):
        """Get the kmer size of the kDataFrameMAP

        :return: kmer size
        :rtype: integer
        """
        pass
