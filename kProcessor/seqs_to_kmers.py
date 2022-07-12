def index(kframe, filename, chunk_size, names_fileName):
    """
    Perform indexing to a sequences file with predefined kmers decoding mode.

    :param kframe: the kDataFrame to be filled with the kmers with their colors
    :type kframe: :class:`kProcessor.kDataFrame`
    :param filename: Sequence(s) file path
    :type filename: str
    :param chunk_size: Number of sequences to parse at once.
    :type chunk_size: int
    :param names_fileName: The TSV names file that contains target sequences headers corresponding to their groups.
    :return: colored kDataFrame with the decoded kmers and their colors.
    :rtype: :class:`kProcessor.colored_kDataFrame`

    Example 1:
        >>> import kProcessor as kp
        >>> KF = kp.kDataFrameMQF(31)
        >>> ckf = kp.index(KF, "seq.fa", 1000, "seq.names")
    
    Example 2:
        >>> import kProcessor as kp
        >>> KF_prot = kp.kDataFramePHMAP(PROTEIN, protein_hasher, {'kSize': 5})
        >>> ckf_prot = kp.index(KF, "seq.fa", 1000, "seq.names")

    """


def countKmersFromFile(kframe, filename, chunk_size):
    """Load the kmers with their counts in the input file into the output kDataframe. Input File can be of formats: fastq,fasta.

    .. note:: The kDataFrame of this function are passed-by-reference. So, it returns nothing.

    :param kframe: the kDataFrame to be filled with the kmers with their counts
    :type kframe: :class:`kProcessor.kDataFrame`
    :param filename: Sequence(s) file path
    :type filename: str
    :param chunk_size: Number of sequences to parse at once.
    :type chunk_size: int

     .. note:: Read more about the usage of parse_params in the FAQ page.

    Example:
        >>> import kProcessor as kp
        >>> KF = kp.kDataFramePHMAP(11)
        >>> kp.parseSequencesFromFile(KF, "seq.fa", 1000) # Fill the KF with the kmers and counts

    """
    pass

def countKmersFromString(seq, kFrame):
    """Load the kmers in the input string into the output kDataframe.

    .. note:: The kDataFrame of this function are passed-by-reference. So, it returns nothing.

    :param kFrame: the kDataFrame to be filled with the kmers with their counts
    :type kFrame: :class:`kProcessor.kDataFrame`
    :param sequence: Sequence to be parsed
    :type sequence: string

    Example:
        >>> import kProcessor as kp
        >>> KF = kDataFramePHMAP(11)
        >>> seq = "ACGATCGATCGATTATATATATCGACGATCGATCGTACGTAGC"
        >>> kp.parseSequencesFromString(seq, KF) # Fill the KF with the kmers and counts

    """
    pass