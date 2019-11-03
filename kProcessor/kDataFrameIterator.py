class kDataFrameIterator:
    """
    Base class for kDataFrame Iterator
    """

    def next(self):
        """
        Increment the iterator to the next kmer
        """

    def getKmer(self):
        """
        Get the kmer at the current iterator position

        :return: Kmer at the current position
        :rtype: string
        """

        pass

    def getHashedKmer(self):
        """
        Get the hash value of the kmer at the current iterator position

        :return: Kmer's hash value at the current position
        :rtype: integer
        """

        pass

    def getCount(self):
        """
        Get the count of the kmer at the current iterator position

        :return: kmer count
        :rtype: integer
        """

        pass

    def setCount(self):
        """
        Sets the count of the current kmer

        :return: True if succeeded, False if failed
        :rtype: boolean
        """

        pass
