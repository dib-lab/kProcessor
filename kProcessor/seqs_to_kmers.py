def index(kframe, parse_params, filename, chunk_size, names_fileName):
    """
    Perform indexing to a sequences file with predefined kmers decoding mode.

    :param kframe: the kDataFrame to be filled with the kmers with their colors
    :type kframe: :class:`kProcessor.kDataFrame`
    :param parse_params: Sequence Decoding parameters
    :type parse_params: dict
    :param filename: Sequence(s) file path
    :type filename: str
    :param chunk_size: Number of sequences to parse at once.
    :type chunk_size: int
    :param names_fileName: The TSV names file that contains target sequences headers corresponding to their groups.
    :return: colored kDataFrame with the decoded kmers and their colors.
    :rtype: :class:`kProcessor.colored_kDataFrame`

     .. note:: Read more about the usage of parse_params in the FAQ page.

    Example:
        >>> import kProcessor as kp
        >>> KF = kp.kDataFrameMQF(31)
        >>> ckf = kp.index(KF,{"mode": 1, "k_size" : 31}, "seq.fa", 1000, "seq.names") # Indexing

    """


def countKmersFromFile(kframe, parse_params, filename, chunk_size):
    """Load the kmers with their counts in the input file into the output kDataframe. Input File can be of formats: fastq,fasta.

    .. note:: The kDataFrame of this function are passed-by-reference. So, it returns nothing.

    :param kframe: the kDataFrame to be filled with the kmers with their counts
    :type kframe: :class:`kProcessor.kDataFrame`
    :param parse_params: Sequence Decoding parameters
    :type parse_params: dict
    :param filename: Sequence(s) file path
    :type filename: str
    :param chunk_size: Number of sequences to parse at once.
    :type chunk_size: int

     .. note:: Read more about the usage of parse_params in the FAQ page.

    Example:
        >>> import kProcessor as kp
        >>> KF = kDataFramePHMAP(11)
        >>> kp.parseSequencesFromFile(KF, "kmers", {"k_size": 11}, "seq.fa", 1000) # Fill the KF with the kmers and counts

    """
    pass

def countKmersFromString(kFrame, parse_params, seq):
    """Load the kmers in the input string into the output kDataframe.

    .. note:: The kDataFrame of this function are passed-by-reference. So, it returns nothing.

    :param kFrame: the kDataFrame to be filled with the kmers with their counts
    :type kFrame: :class:`kProcessor.kDataFrame`
    :param parse_params: Sequence Decoding parameters
    :type parse_params: dict
    :param sequence: Sequence to be parsed
    :type sequence: string

     .. note:: Read more about the usage of parse_params in the FAQ page.

    Example:
        >>> import kProcessor as kp
        >>> KF = kDataFramePHMAP(11)
        >>> seq = "ACGATCGATCGATTATATATATCGACGATCGATCGTACGTAGC"
        >>> kp.parseSequencesFromString(KF, "kmers", {"k_size": 11}, seq) # Fill the KF with the kmers and counts

    """
    pass