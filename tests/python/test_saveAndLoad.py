import os
import shutil
import unittest
import random
import kProcessor as kp
from params import test_params


class TestSaveAndLoad(unittest.TestCase):
    ksize = 21

    params_PHMAP = test_params(ksize, "PHMAP")
    params_MAP = test_params(ksize, "MAP")

    def test_MQF(self):

        my_tmpdir = "tmp" + str(random.randint(1, 9999))
        os.mkdir(my_tmpdir)

        params = test_params(self.ksize, "MQF")
        params.Hasher = kp.IntegerHasher(self.ksize)

        KF = params.new_kf(self.ksize)
        self.assertTrue(KF.empty())

        # Create random kmers list
        kmers_list = params.generate_kmers(kmers_no=20)
        kmers_hashes = [params.Hasher.hash(kmer[0]) for kmer in kmers_list]

        for kmer in kmers_list:
            KF.insert(kmer[0], kmer[1])

        fileName = os.path.join(my_tmpdir, "tmp.kdataframeMQF")
        KF.save(fileName)
        KF_LOADED = kp.kDataFrame.load(fileName)

        it = KF_LOADED.begin()
        count = 0

        while it != KF_LOADED.end():
            self.assertTrue(it.getHashedKmer() in kmers_hashes)
            count += 1
            it.next()

        self.assertEqual(count, len(kmers_list))
        shutil.rmtree(my_tmpdir)

    def test_MAP(self):

        my_tmpdir = "tmp" + str(random.randint(1, 9999))
        os.mkdir(my_tmpdir)

        params = test_params(self.ksize, "MAP")
        params.Hasher = kp.TwoBitsHasher(self.ksize)

        KF = params.new_kf(self.ksize)
        self.assertTrue(KF.empty())

        # Create random kmers list
        kmers_list = params.generate_kmers(kmers_no=20)
        kmers_hashes = [params.Hasher.hash(kmer[0]) for kmer in kmers_list]

        for kmer in kmers_list:
            KF.insert(kmer[0], kmer[1])

        fileName = os.path.join(my_tmpdir, "tmp.kdataframeMQF")
        KF.save(fileName)
        KF_LOADED = kp.kDataFrame.load(fileName)

        it = KF_LOADED.begin()
        count = 0

        while it != KF_LOADED.end():
            self.assertTrue(it.getHashedKmer() in kmers_hashes)
            count += 1
            it.next()

        self.assertEqual(count, len(kmers_list))
        shutil.rmtree(my_tmpdir)

    def test_PHMAP(self):

        my_tmpdir = "tmp" + str(random.randint(1, 9999))
        os.mkdir(my_tmpdir)

        params = test_params(self.ksize, "PHMAP")
        params.Hasher = kp.TwoBitsHasher(self.ksize)

        KF = params.new_kf(self.ksize)
        self.assertTrue(KF.empty())

        # Create random kmers list
        kmers_list = params.generate_kmers(kmers_no=20)
        kmers_hashes = [params.Hasher.hash(kmer[0]) for kmer in kmers_list]

        for kmer in kmers_list:
            KF.insert(kmer[0], kmer[1])

        fileName = os.path.join(my_tmpdir, "tmp.kdataframeMQF")
        KF.save(fileName)
        KF_LOADED = kp.kDataFrame.load(fileName)

        it = KF_LOADED.begin()
        count = 0

        while it != KF_LOADED.end():
            self.assertTrue(it.getHashedKmer() in kmers_hashes)
            count += 1
            it.next()

        self.assertEqual(count, len(kmers_list))
        shutil.rmtree(my_tmpdir)
