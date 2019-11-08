import os
import shutil
import unittest
import random
import kProcessor as kp
from params import test_params


class TestKDataFramePHMAP(unittest.TestCase):
    params = test_params(kSize=21, kDataFrameType="PHMAP")
    params.Hasher = kp.TwoBitsHasher(21)

    def test_emptykDataFrame(self):
        framesToBeTested = []
        kSizes = [21, 31]
        for kSize in kSizes:
            framesToBeTested.append(self.params.new_kf(kSize))

        for kFrame in framesToBeTested:
            self.assertTrue(kFrame.empty())
            self.assertEqual(kFrame.size(), 0)

    def test_insertOneTime(self):
        print(self._testMethodName)
        kSize = 31
        _kmer, rand_count = self.params.generate_singleKmer()

        kf = self.params.create_new_kf()
        self.assertTrue(kf.empty())

        kf.insert(_kmer, rand_count)
        self.assertFalse(kf.empty())
        self.assertEqual(kf.getCount(_kmer), rand_count)

    def test_insertNTimes(self):
        print(self._testMethodName)

        # Create random kmers list
        kmers_list = self.params.generate_kmers(kmers_no=20)

        # Empty kDataFrames
        kFrame = self.params.create_new_kf()
        self.assertTrue(kFrame.empty())

        insertedKmers = 0
        kmers_hash_values = [self.params.Hasher.hash(kmer[0]) for kmer in kmers_list]

        # Insert all kmers
        for kmer in kmers_list:
            if kFrame.insert(kmer[0], kmer[1]):
                insertedKmers += 1

        # Assert all kmers are inserted
        self.assertEqual(insertedKmers, len(kmers_list))

        # Verify inserted kmers
        it = kFrame.begin()
        while it != kFrame.end():
            self.assertTrue(it.getHashedKmer() in kmers_hash_values)
            it.next()

        # Verify all inserted kmers
        for kmer in kmers_list:
            c = kFrame.getCount(kmer[0])
            self.assertEqual(c, kmer[1])

    def test_eraseKmers(self):
        print(self._testMethodName)

        # Create random kmers list
        kmers_list = self.params.generate_kmers(kmers_no=20)

        # Empty kDataFrames
        kFrames = self.params.create_empty_kframes(2)
        self.assertTrue(kFrames[0].empty())
        self.assertTrue(kFrames[1].empty())

        # Insert kmers
        for i in range(len(kFrames)):
            for kmer in kmers_list:
                self.assertTrue(kFrames[i].insert(kmer[0], kmer[1]))

        # Erasing all kmers
        for i in range(len(kFrames)):
            for kmer in kmers_list:
                self.assertTrue(kFrames[i].erase(kmer[0]))

        # Check that all kmers have been erased
        for i in range(len(kFrames)):
            for kmer in kmers_list:
                self.assertEqual(kFrames[i].getCount(kmer[0]), 0)

    def test_iterateOverAllKmers(self):
        print(self._testMethodName)

        # Create random kmers list
        kmers_list = self.params.generate_kmers(kmers_no=20)

        # Empty kDataFrames
        kFrames = self.params.create_empty_kframes(2)
        self.assertTrue(kFrames[0].empty())
        self.assertTrue(kFrames[1].empty())

        # Insert kmers
        for i in range(len(kFrames)):
            for kmer in kmers_list:
                self.assertTrue(kFrames[i].insert(kmer[0], kmer[1]))

        # Get all inserted kmers counts
        inserted_counts = set([kmer[1] for kmer in kmers_list])

        for kFrame in kFrames:
            it = kFrame.begin()
            kframe_kmers_counts = set()
            while it != kFrame.end():
                count = it.getCount()
                kframe_kmers_counts.add(count)
                it.next()

            self.assertEqual(len(kframe_kmers_counts.intersection(inserted_counts)), len(inserted_counts))

    def test_saveAndIterateOverAllKmers(self):
        my_tmpdir = "tmp" + str(random.randint(1, 9999))
        os.mkdir(my_tmpdir)
        print(self._testMethodName)

        # Create random kmers list
        kmers_list = self.params.generate_kmers(kmers_no=20)

        # Empty kDataFrames
        kFrames = self.params.create_empty_kframes(2)
        self.assertTrue(kFrames[0].empty())
        self.assertTrue(kFrames[1].empty())

        # Insert kmers
        for i in range(len(kFrames)):
            for kmer in kmers_list:
                self.assertTrue(kFrames[i].insert(kmer[0], kmer[1]))

        # Get all inserted kmers counts
        inserted_counts = set([kmer[1] for kmer in kmers_list])
        kmers_hash_values = [self.params.Hasher.hash(kmer[0]) for kmer in kmers_list]

        fileName = os.path.join(my_tmpdir, "tmp.kdataframe")

        for i in range(len(kFrames)):
            kFrames[i].save(fileName + "_" + str(i))

        loaded_kFrames = []
        for i in range(len(kFrames)):
            loaded_kFrames.append(kp.kDataFrame.load(fileName + "_" + str(i)))

        for kFrame in loaded_kFrames:
            it = kFrame.begin()
            kframe_kmers_counts = set()
            while (it != kFrame.end()):
                count = it.getCount()
                kframe_kmers_counts.add(count)
                self.assertTrue(it.getHashedKmer() in kmers_hash_values)
                self.assertTrue(it.getCount() in inserted_counts)
                it.next()

        shutil.rmtree(my_tmpdir)
