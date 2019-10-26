def kFrameUnion(input):
    """
    Calculate the union of the kDataFrames. The result kDataframe will have all the kmers in the input list of kDataframes.
    The count of the kmers equals to the sum of the kmer count in the input list.

    .. warning:: This function works only with :class:`kProcessor.kDataFrameMQF`.

    :param input: List of kDataFrames
    :type input: list of :class:`kProcessor.kDataFrameMQF`
    :return: New kDataFrame object holding the union of kmers in the kDataFrames list.
    :rtype: :class:`kProcessor.kDataFrame`

    """

def kFrameIntersect(input):
    """
    Calculate the intersect of the kDataFrames. The result kDataframe will have only kmers that exists in all the kDataframes.
    The count of the kmers equals to the minimum of the kmer count in the input list.

    .. warning:: This function works only with :class:`kProcessor.kDataFrameMQF`.

    :param input: List of kDataFrames
    :type input: list of :class:`kProcessor.kDataFrameMQF`
    :return: New kDataFrame object holding the intersection of kmers in the kDataFrames list.
    :rtype: :class:`kDataFrame`

    """

def kFrameDiff(input):
    """
    Calculate the difference of the kDataframes.
    The result kDataframe will have only kmers that exists in the first kDataframe and not in any of the rest input kDataframes.
    The count of the kmers equals to the count in the first kDataframe.

    .. warning:: This function works only with :class:`kProcessor.kDataFrameMQF`.

    :param input: List of kDataFrames
    :type input: list of :class:`kProcessor.kDataFrameMQF`
    :return: New kDataFrame object holding the difference of kmers in the kDataFrames list.
    :rtype: :class:`kDataFrame`

    """
