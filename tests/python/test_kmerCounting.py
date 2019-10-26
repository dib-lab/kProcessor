import os
import shutil
import unittest
import random
import kProcessor as kp
from params import test_params


class TestkmerCounting(unittest.TestCase):
    ksize = 21

    def test_kmerCountingFASTQ(self):
        params_MQF = test_params(self.ksize, "MQF")
        params_PHMAP = test_params(self.ksize, "PHMAP")
        params_MAP = test_params(self.ksize, "MAP")

        KF_MQF = params_MQF.new_kf(self.ksize)
        KF_PHMAP = params_PHMAP.new_kf(self.ksize)
        KF_MAP = params_MAP.new_kf(self.ksize)

        # kmerCounting
        kp.parseSequencesFromFile(KF_MQF, "kmers", {"k_size": self.ksize}, params_MQF.fastqFiles[0], 1000)
        kp.parseSequencesFromFile(KF_PHMAP, "kmers", {"k_size": self.ksize}, params_PHMAP.fastqFiles[0], 1000)
        kp.parseSequencesFromFile(KF_MAP, "kmers", {"k_size": self.ksize}, params_MAP.fastqFiles[0], 1000)

        self.assertTrue(KF_MQF.size())
        self.assertTrue(KF_PHMAP.size())
        self.assertTrue(KF_MAP.size())

    def test_kmerCountingFASTA(self):
        params_MQF = test_params(self.ksize, "MQF")
        params_PHMAP = test_params(self.ksize, "PHMAP")
        params_MAP = test_params(self.ksize, "MAP")

        KF_MQF = params_MQF.new_kf(self.ksize)
        KF_PHMAP = params_PHMAP.new_kf(self.ksize)
        KF_MAP = params_MAP.new_kf(self.ksize)

        self.assertFalse(KF_MQF.size())
        self.assertFalse(KF_PHMAP.size())
        self.assertFalse(KF_MAP.size())

        # kmerCounting
        kp.parseSequencesFromFile(KF_MQF, "kmers", {"k_size": self.ksize}, params_MQF.fasta_file, 1000)
        kp.parseSequencesFromFile(KF_PHMAP, "kmers", {"k_size": self.ksize}, params_PHMAP.fasta_file, 1000)
        kp.parseSequencesFromFile(KF_MAP, "kmers", {"k_size": self.ksize}, params_MAP.fasta_file, 1000)

        self.assertTrue(KF_MQF.size())
        self.assertTrue(KF_PHMAP.size())
        self.assertTrue(KF_MAP.size())
