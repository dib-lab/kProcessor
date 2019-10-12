from kProcessor.kDataFrame import kDataFrame

class colored_kDataFrame(kDataFrame):
    """colored_kDataFrame class

    .. note:: the colored_kDataFrame Inherits all the functions from :class:`kProcessor.kDataFrame` plus other new functions.

    *Introduction*:
        - The colored_kDataFrame class holds the Kmers colors instead of their count.
        - The **color** is an integer represents the targets which contains that kmer.

        Example:
            **color:** ``1``: represents the transcripts ``transcript_A`` , ``transcript_B`` and ``transcript_C``
            **color:** ``2``: represents the transcripts ``transcript_A`` , ``transcript_B``

            **kmer:** ``ACTGATCGATCGTACGAC`` has the **color** `2`, that means it's found in both `transcript_A` and `transcript_B`
            **kmer:** ``ATAAGCATTTACAGCAAT`` has the **color** `1`, that means it's found in both `transcript_A` , `transcript_B` and `transcript_C`


    """

    pass


    def getKmerColor(self, kmer):
        """
        Get the color of the kmer

        :param kmer: Kmer string
        :type kmer: str
        :return: The color of the kmer
        :rtype: int
        """
        pass


    def getSamplesIDForKmer(self, kmer):
        """
        Get all sample IDs that contains that kmer.


        :param kmer: Kmer string
        :type kmer: str
        :return: List of all samples IDs associated with that kmer.
        :rtype: list
        """

    def getSamplesIDForColor(self, color):
        """
        Get all sample IDs that contains that kmer.


        :param color: Kmer color
        :type color: int
        :return: List of all samples IDs associated with that color.
        :rtype: list
        """


    def names_map(self):
        """
        Get the names map dictionary that represents sample ID as key and its group name as value.

        :return: names map dictionary.
        :rtype: dict
        """

    def inverse_names_map(self):
        """
        Get the names map dictionary that represents group name as key and its sample ID as value.

        :return: inverse names map dictionary.
        :rtype: dict
        """

