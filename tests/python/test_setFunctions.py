import os
import shutil
import unittest
import random
import kProcessor as kp
from params import test_params


class TestKDataFrameMQF(unittest.TestCase):
    params = test_params(kSize=21, kDataFrameType="MQF")
    params.Hasher = kp.IntegerHasher(21)

    def test_kFrameUnion(self):
        print(self._testMethodName)

        # Empty kDataFrames
        kFrames_vec = self.params.create_empty_kframes(2)
        self.assertTrue(kFrames_vec[0].empty())
        self.assertTrue(kFrames_vec[1].empty())

        # Create random kmers list
        kmers_list1 = self.params.generate_kmers(kmers_no=20)
        inserted_counts1 = set([kmer[1] for kmer in kmers_list1])

        kmers_list2 = self.params.generate_kmers(kmers_no=20)
        inserted_counts2 = set([kmer[1] for kmer in kmers_list2])

        # Inserting Kmers
        inserted_kmers_hashes_1 = [self.params.Hasher.hash(kmer[0]) for kmer in kmers_list1]
        inserted_kmers_hashes_2 = [self.params.Hasher.hash(kmer[0]) for kmer in kmers_list2]
        inserted_kmers_hashes_1 = set(inserted_kmers_hashes_1)
        inserted_kmers_hashes_2 = set(inserted_kmers_hashes_2)
        union_kmers_hashes = inserted_kmers_hashes_1.union(inserted_kmers_hashes_2)

        for kmer in kmers_list1:
            self.assertTrue(kFrames_vec[0].insert(kmer[0], kmer[1]))

        for kmer in kmers_list2:
            self.assertTrue(kFrames_vec[1].insert(kmer[0], kmer[1]))

        # Apply kFrameUnion
        union_kFrame = kp.kFrameUnion(kFrames_vec)

        # Total kmers extracted from union kFrames_vec
        kmers_count = 0

        it = union_kFrame.begin()
        while it != union_kFrame.end():
            self.assertTrue(it.getHashedKmer() in union_kmers_hashes)
            kmers_count += 1
            it.next()

        self.assertEqual(kmers_count, len(union_kmers_hashes))

    def test_kFrameIntersect(self):
        print(self._testMethodName)

        # Empty kDataFrames
        kFrames_vec = self.params.create_empty_kframes(2)
        self.assertTrue(kFrames_vec[0].empty())
        self.assertTrue(kFrames_vec[1].empty())

        # Create random kmers list
        kmers_list1 = self.params.generate_kmers(kmers_no=20)
        inserted_counts1 = set([kmer[1] for kmer in kmers_list1])

        kmers_list2 = self.params.generate_kmers(kmers_no=20)
        # Replicate some kmers from kmers_list1 in kmers_list2 to make sure len(intersection) > 0
        kmers_list2 += kmers_list1[0:10]

        inserted_counts2 = set([kmer[1] for kmer in kmers_list2])

        # Inserting Kmers
        inserted_kmers_hashes_1 = [self.params.Hasher.hash(kmer[0]) for kmer in kmers_list1]
        inserted_kmers_hashes_2 = [self.params.Hasher.hash(kmer[0]) for kmer in kmers_list2]
        inserted_kmers_hashes_1 = set(inserted_kmers_hashes_1)
        inserted_kmers_hashes_2 = set(inserted_kmers_hashes_2)
        intersect_kmers_hashes = inserted_kmers_hashes_1.intersection(inserted_kmers_hashes_2)

        for kmer in kmers_list1:
            self.assertTrue(kFrames_vec[0].insert(kmer[0], kmer[1]))

        for kmer in kmers_list2:
            self.assertTrue(kFrames_vec[1].insert(kmer[0], kmer[1]))

        # Apply kFrameUnion
        intersect_kFrame = kp.kFrameIntersect(kFrames_vec)

        # Total kmers extracted from union kFrames_vec
        kmers_count = 0

        it = intersect_kFrame.begin()
        while it != intersect_kFrame.end():
            self.assertTrue(it.getHashedKmer() in intersect_kmers_hashes)
            kmers_count += 1
            it.next()


        self.assertEqual(kmers_count, len(intersect_kmers_hashes))

    def test_kFrameDiff(self):
        print(self._testMethodName)

        # Empty kDataFrames
        kFrames_vec = self.params.create_empty_kframes(2)
        self.assertTrue(kFrames_vec[0].empty())
        self.assertTrue(kFrames_vec[1].empty())

        # Create random kmers list
        kmers_list1 = self.params.generate_kmers(kmers_no=20)
        inserted_counts1 = set([kmer[1] for kmer in kmers_list1])

        kmers_list2 = self.params.generate_kmers(kmers_no=20)
        # Replicate some kmers from kmers_list1 in kmers_list2 to make sure len(intersection) > 0
        kmers_list2 += kmers_list1[0:10]

        inserted_counts2 = set([kmer[1] for kmer in kmers_list2])

        # Inserting Kmers
        inserted_kmers_hashes_1 = [self.params.Hasher.hash(kmer[0]) for kmer in kmers_list1]
        inserted_kmers_hashes_2 = [self.params.Hasher.hash(kmer[0]) for kmer in kmers_list2]
        inserted_kmers_hashes_1 = set(inserted_kmers_hashes_1)
        inserted_kmers_hashes_2 = set(inserted_kmers_hashes_2)
        difference_kmers_hashes = inserted_kmers_hashes_1.difference(inserted_kmers_hashes_2)

        for kmer in kmers_list1:
            self.assertTrue(kFrames_vec[0].insert(kmer[0], kmer[1]))

        for kmer in kmers_list2:
            self.assertTrue(kFrames_vec[1].insert(kmer[0], kmer[1]))

        # Apply kFrameUnion
        diff_kFrame = kp.kFrameDiff(kFrames_vec)

        # Total kmers extracted from union kFrames_vec
        kmers_count = 0

        it = diff_kFrame.begin()
        while it != diff_kFrame.end():
            self.assertTrue(it.getHashedKmer() in difference_kmers_hashes)
            kmers_count += 1
            it.next()

        self.assertEqual(kmers_count, len(difference_kmers_hashes))
