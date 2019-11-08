import os
import shutil
import unittest
import random
import kProcessor as kp
from params import test_params
import pytest


class TestkmerDecoder(unittest.TestCase):
    generic_params = test_params(10, "PHMAP")
    seq = "ACGTAGCATGCATGACGATGCTAGCGTGATGCTAGCTAGTCAGTAGCATGC"

    def test_kmers(self):
        KD_params = {
            "mode" : 1,
            "k_size": 10,
        }


        KF_seq = kp.kDataFramePHMAP(KD_params["k_size"])
        KF_file = kp.kDataFramePHMAP(KD_params["k_size"])


        # void countKmersFromFile(kDataFrame * output, string mode, std::map<std::string, int> params, string filename, int chunk_size = 1000);
        kp.countKmersFromFile(KF_file, KD_params, self.generic_params.small_fasta_file, 1)
        self.assertFalse(KF_file.empty())

        kp.countKmersFromString(KF_seq, KD_params, self.seq)
        self.assertFalse(KF_seq.empty())


    def test_skipmers(self):
        KD_params = {
            "mode": 2,
            "k_size": 10,
            "m": 2,
            "n": 3
        }


        KF_seq = kp.kDataFramePHMAP(KD_params["k_size"])
        KF_file = kp.kDataFramePHMAP(KD_params["k_size"])

        kp.countKmersFromFile(KF_file, KD_params, self.generic_params.small_fasta_file, 1)
        self.assertFalse(KF_file.empty())

        kp.countKmersFromString(KF_seq, KD_params, self.seq)
        self.assertFalse(KF_seq.empty())

        self.assertEqual(KF_seq.size(), KF_file.size())

    # def test_minimizers(self):
    #     KD_params = {
    #         "mode": "minimizers",
    #         "params": {
    #             "k_size": 5,
    #             "w": 10,
    #         }
    #     }
    #
    #     KD_seq = kp.initialize_kmerDecoder(KD_params["mode"], KD_params["params"])
    #     KF_seq = kp.kDataFramePHMAP(KD_params["params"]["k_size"])
    #
    #
    #     KD_file = kp.initialize_kmerDecoder(self.generic_params.small_fasta_file, 1, KD_params["mode"], KD_params["params"])
    #     KF_file = kp.kDataFramePHMAP(KD_params["params"]["k_size"])
    #
    #     kp.parseSequences(KD_file, KF_file)
    #     self.assertFalse(KF_file.empty())
    #
    #     kp.countKmersFromString(KD_seq, self.seq, KF_seq)
    #     self.assertFalse(KF_seq.empty())
    #
    #     self.assertEqual(KF_seq.size(), KF_file.size())


