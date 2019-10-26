def index(kFrame, mode, params, filename, chunk_size, names_fileName):
    """
    Perform indexing to a sequences file with predefined kmers decoding mode.

    :param kFrame: the kDataFrame to be filled with the kmers with their colors
    :type kFrame: :class:`kProcessor.kDataFrame`
    :param mode: String defines the mode of kmerDecoder to parse sequences.
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
    """


def parseSequencesFromFile(kFrame, mode, params, filename, chunk_size):
    """Load the kmers in the input file into the output kDataframe. Input File can be of formats: fastq,fasta.

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

    """
    pass

#void parseSequencesFromString(kDataFrame * frame, string mode, std::map<std::string, int> params, string sequence);

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

    """
    pass

def initialize_kmerDecoder(*args):
    """
    Initialize a kmerDecoder object

    .. note:: This function has multiple overridden functions.

    kmerDecoder is important to tell the kProcessor how to parse sequences.
    It can be initialized to parse sequences files "FASTA/Q" or sequence strings.


    :param filename: File path
    :type filename: str

    :param chunkSize: Number of sequences to parse at once.
    :type chunkSize: int

    :param mode: String defines the mode of kmerDecoder to parse sequences.
    :type mode: str

    :param params: kmerDecoder parameters
    :type params: dict

    :param kSize: Kmer size
    :type kSize: int

    :param hash_mode: The hash function will be used by the kmerDecoder for kmers hashing.
    :type hash_mode: int

    :return: kmerDecoder Object
    :rtype: kmerDecoder

    ---

    *Modes*
        There are three modes for decoding the sequences substrings ``kmers``, ``skipmers``, ``minimizers``

    *Modes Parameters*
        Mode parameters is a dictionary contains the parameters to initialize the kmerDecoder instance.

        1. ``kmers`` mode parameters : {"k_size" : ``k`` } , where k is an integer.
        2. ``skipmers`` mode parameters : {"k_size" : ``k``, "m" : ``m``, "n" : ``n`` } , where k,m,n are integers.
        3. ``minimizers`` mode parameters : {"k_size" : ``k``, "w" : ``w`` } , where k,w are integers.

    ---

    Example 1:
        *Initializing kmerDecoder to parse sequences from FASTA/FASTQ file using the popular mode [kmers] with kmerSize 31*

        *Function Signature*:
            `initialize_kmerDecoder(filename, chunkSize, mode, params)`

        >>> import kProcessor as kp
        >>> file_path = "path/to/seq.fa"
        >>> mode = "kmers"
        >>> kSize = 31
        >>> params = {"k_size" : kSize}
        >>> chunk_size = 1000 # Load and parse 1000 sequences per time
        >>> KD = kp.initialize_kmerDecoder(file_path, chunk_size, mode, params)
        >>> KF = kDataFrameMQF(kSize) # The kmer size of the kDataFrame and kmerDecoder MUST be the same.
        >>> kp.parseSequences(kD, KF) # Now the KF has been filled with the kmers and their counts.

    Example 2:
        *Initialize kmerDecoder to parse a DNA string using the skipmers mode*

        *Function Signature*:
            `initialize_kmerDecoder(mode, params)`


        >>> import kProcessor as kp
        >>> seq = "ACGTACGTACGATCGTAGCATCGATTACGAATCGATCGATA"
        >>> mode = "skipmers"
        >>> k, m, n = 10, 2, 3
        >>> params = {"k_size" : k, "m" : m, "n" : n}
        >>> KD = kp.initialize_kmerDecoder(mode, params)
        >>> KF = kDataFramePHMAP(10)
        >>> kp.parseSequencesFromString(KD, seq, KF)

    Example 3:
        *Initialize kmerDecoder without parsing any sequences, to be used in other functions*

        *Hashing Modes*:
             Mode 0: Murmar Hashing  | Irreversible
             Mode 1: Integer Hashing | Reversible
             Mode 2: TwoBitsHashing  | Not considered hashing, just the two bits representation of the kmer

        *Function Signature*
            `initialize_kmerDecoder(kSize, hash_mode)`

        >>> import kProcessor as kp
        >>> KD = initialize_kmerDecoder(31, 1) # The most common used mode is the Integer Hashing

    """
    pass


def kmerDecoder_setHashing(KD, hash_mode, canonical):
    """
    set the hashing mode of kmerDecoder object

    .. note:: the ``KD`` parameter is passed by reference.

    .. note:: The main usage of this function is to set KD to hash kmers in non-canonical mode.

    :param KD: kmerDecoder object to set its hashing mode
    :type KD: kmerDecoder
    :param hash_mode: The hashing mode to be used in kmers hashing
    :type hash_mode: int
    :param canonical: Determining whether it should hash the canonical representation of the kmers or not.
    :type canonical: bool

    *Hashing Modes*:
             Mode 0: Murmar Hashing  | Irreversible
             Mode 1: Integer Hashing | Reversible
             Mode 2: TwoBitsHashing  | Not considered hashing, just the two bits representation of the kmer

    **Example:**

        >>> import kProcessor as kp
        >>> KD = kp.initialize_kmerDecoder(31, 1)
        >>> kp.kmerDecoder_setHashing(KD, 2)

    """