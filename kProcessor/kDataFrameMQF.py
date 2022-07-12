
from kProcessor.kDataFrame import kDataFrame

class kDataFrameMQF(kDataFrame):
    """The abstract base class defining a kDataFrameMQF.
    """

    def __init__(self, kSize):
        """Instantiate a kDataFrameMQF object with predefined kmer size.

        :param kSize: Kmer Size
        :type kSize: integer
        :return: :class:`kProcessor.kDataFrameMQF`

        Instantiation Example:
            >>> import kProcessor as kp
            >>> KF_MQF_1 = kp.kDataFrameMQF(31) # kSize = 31
            >>> KF_MQF_1 = kp.kDataFrameMQF(SKIPMERS, integer_hasher, {'m': 2, 'n': 3, 'k': 10}) # Reading mode = skipmers, hashing mode = integer hashing, (m, n, k) are the skipmers params.
            >>> KF_MQF_2 = kp.kDataFrameMQF(PROTEIN, protein_hasher, {'kSize': 5}); # Reading/hashing mode = protein, kSize = 5
            >>> KF_MQF_3 = kp.kDataFrameMQF(PROTEIN, proteinDayhoff_hasher, {'kSize': 11}); # Reading mode = protein, hashing mode = dayhoff encoding, kSize = 11

            
        .. note:: Read more about hashing modes in the FAQ page.

        """
        pass


    def getTwin(self):
        """creates a new ``kDataFrameMQF`` using the same parameters as the current ``kDataFrameMQF``.
    
        :return: A shallow copy of the current ``kDataFrameMQF``.
        :rtype: kDataFrameMQF
    
        """

    def reserve(self, n):
        """Request a capacity change so that the kDataFrameMQF can approximately hold at least n kmers

        :param n: Minimum number of kmers
        :type n: integer

         .. note:: Read more about the usage of `reserve(n)` in the FAQ page.

        """
        pass

    def insert(self, kmer, N = 1):
        """Insert the kmer N time in the kDataFrameMQF, or increment the kmer count with N if it is already exists.

        :param kmer: The Kmer to increment its count
        :type kmer: string
        :param N: Kmer count (Optional, Default = 1)
        :type N: integer
        :return: Boolean value indicating whether the kmer is inserted or not
        :rtype: bool
        """

        pass

    def setCount(self, kmer, N):
        """Set the kmer's count to N time in the kDataFrameMQF

        :param kmer: The Kmer to set its count
        :type kmer: string
        :param N: Kmer count
        :type N: integer
        :return: Boolean value indicating whether the kmer is inserted or not
        :rtype: bool
        """
        pass

    def getCount(self, kmer):
        """Retrieve number of times the kmer was inserted in the kDataFrameMQF

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
        :rtype: bool
        """
        pass

    def size(self):
        """Number of kmers in the kDataFrameMQF

        :return: The number of kmers in the kDataFrameMQF
        :rtype: integer
        """
        pass

    def max_size(self):
        """Maximum number of kmers that the kDataFrameMQF can hold.

        :return: The maximum number of kmers that the kDataFrameMQF can hold.
        :rtype: integer
        """

        pass

    def empty(self):
        """Check whether the kDataFrameMQF is empty of kmers or not.

        :return: Boolean value indicating whether the kDataFrameMQF is empty, i.e. whether its size is 0
        :rtype: boolean

        """
        pass

    def load_factor(self):
        """Retrieving the current load factor of the kDataFrameMQF in percentage to indicate how full is it.

        :return: The current load factor in the kDataFrameMQF.
        :rtype: integer
        """
        pass

    def max_load_factor(self):
        """ Retrieving the maximum load factor of the kDataFrameMQF in percentage.

        :return: The maximum load factor in the kDataFrameMQF.
        :rtype: integer
        """
        pass

    def begin(self):
        """ Instantiate a kDataFrameIterator object pointing to the first kmer position
        :return: An iterator at the begin of the kDataFrameMQF.
        :rtype: :class:`kProcessor.kDataFrameIterator`
        """

        pass

    def end(self):
        """ Instantiate a kDataFrameIterator object pointing to the last kmer position

        :return: An iterator at the end of the kDataFrameMQF.
        :rtype: :class:`kProcessor.kDataFrameIterator`

        """
        pass

    def save(self):
        """
        Serialize the kDataFrameMQF on the disk in a form of binary file alongside other metadata files.
        """

        pass


    @staticmethod
    def load(filePath):
        """ A static method to load a kDataFrameMQF file from disk.

        .. note:: Load the file without the extension [.mqf, .map, .phmap]

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
        """Get the kmer size of the kDataFrameMQF

        :return: kmer size
        :rtype: integer
        """
        pass
