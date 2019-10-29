def index(kFrame, mode, params, filename, chunk_size, names_fileName):
    """
    Perform indexing to a sequences file with predefined kmers decoding mode.

    :param kFrame: the kDataFrame to be filled with the kmers with their colors
    :type kFrame: :class:`kProcessor.kDataFrame`
    :param mode: String defines the mode of to parse sequences.
    :type mode: str
    :param params: Sequence Decoding parameters
    :type params: dict
    :param filename: Sequence(s) file path
    :type filename: str
    :param chunk_size: Number of sequences to parse at once.
    :type chunk_size: int
    :param names_fileName: The TSV names file that contains target sequences headers corresponding to their groups.
    :return: colored kDataFrame with the decoded kmers and their colors.
    :rtype: :class:`kProcessor.colored_kDataFrame`

    Example:
        >>> import kProcessor as kp
        >>> KF = kp.kDataFrameMQF(31)
        >>> ckf = kp.index(KF, "kmers", {"k_size" : 31}, "seq.fa", 1000, "seq.names") # Indexing

    """


def parseSequencesFromFile(kFrame, mode, params, filename, chunk_size):
    """Load the kmers with their counts in the input file into the output kDataframe. Input File can be of formats: fastq,fasta.

    .. note:: The kDataFrame of this function are passed-by-reference. So, it returns nothing.

    :param kFrame: the kDataFrame to be filled with the kmers with their counts
    :type kFrame: :class:`kProcessor.kDataFrame`
    :param mode: String defines the mode of kmerDecoder to parse sequences.
    :type mode: str
    :param params: Sequence Decoding parameters
    :type params: dict
    :param filename: Sequence(s) file path
    :type filename: str
    :param chunk_size: Number of sequences to parse at once.
    :type chunk_size: int

    Example:
        >>> import kProcessor as kp
        >>> KF = kDataFramePHMAP(11)
        >>> kp.parseSequencesFromFile(KF, "kmers", {"k_size": 11}, "seq.fa", 1000) # Fill the KF with the kmers and counts

    """
    pass

def parseSequencesFromString(kFrame, mode, params, seq):
    """Load the kmers in the input string into the output kDataframe.

    .. note:: The kDataFrame of this function are passed-by-reference. So, it returns nothing.

    :param kFrame: the kDataFrame to be filled with the kmers with their counts
    :type kFrame: :class:`kProcessor.kDataFrame`
    :param mode: String defines the mode of kmerDecoder to parse sequences.
    :type mode: str
    :param params: Sequence Decoding parameters
    :type params: dict
    :param sequence: Sequence to be parsed
    :type sequence: string

    Example:
        >>> import kProcessor as kp
        >>> KF = kDataFramePHMAP(11)
        >>> seq = "ACGATCGATCGATTATATATATCGACGATCGATCGTACGTAGC"
        >>> kp.parseSequencesFromString(KF, "kmers", {"k_size": 11}, seq) # Fill the KF with the kmers and counts


    """
    pass