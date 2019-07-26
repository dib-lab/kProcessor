import random
import os
import kProcessor as kp


class test_params:
    kSize = 31
    kDataFrameType = str()
    kDataFrames = list()
    new_kf = None
    Hasher = None

    def __init__(self, kSize, kDataFrameType):
        self.kDataFrameType = kDataFrameType
        self.kSize = kSize

        # Set Absolute Path
        self.abs_path = os.path.normpath(os.getcwd() + os.sep + os.pardir)
        self.test_data_path = os.path.join(self.abs_path, "testData")

        # Fasta file and Names file
        self.fasta_file = os.path.join(self.test_data_path, "test1.fa")
        self.names_file = os.path.join(self.test_data_path, "test1.fa.names")

        fastq_files = ["test.noN.fastq", "test2.noN.fastq", "test2.noN.fastq"]
        self.fastqFiles = [os.path.join(self.test_data_path, fastq) for fastq in fastq_files]

        if self.kDataFrameType == "MQF":
            self.new_kf = kp.kDataFrameMQF
        elif self.kDataFrameType == "MAP":
            self.new_kf = kp.kDataFrameMAP
        elif self.kDataFrameType == "PHMAP":
            self.new_kf = kp.kDataFramePHMAP

        # self.test_build_kDataFrames()

    def create_new_kf(self):
        return self.new_kf(self.kSize)

    def test_build_kDataFrames(self):

        for fastq in self.fastqFiles:
            KD = kp.initialize_kmerDecoder(fastq, 1000, "kmers", {"k_size": self.kSize})
            KF = self.new_kf(self.kSize)
            kp.parseSequences(KD, KF)
            self.kDataFrames.append(KF)

    def create_empty_kframes(self, N):
        kframes = []
        for i in range(N):
            kframes.append(self.create_new_kf())

        return kframes

    def generate_singleKmer(self):
        _kmer = "".join([random.choice('ACGT') for _ in range(self.kSize)])
        _count = random.randint(1, 1000)
        return _kmer, _count

    def generate_kmers(self, kmers_no):
        kmers_list = []
        for i in range(kmers_no):
            _kmer = "".join([random.choice('ACGT') for _ in range(self.kSize)])
            _count = random.randint(1, 1000)
            kmers_list.append((_kmer, _count))

        return kmers_list
