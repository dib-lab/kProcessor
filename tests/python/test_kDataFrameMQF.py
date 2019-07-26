import os
import shutil
import unittest
import random
import kProcessor as kp
from params import test_params


class TestKDataFrameMQF(unittest.TestCase):
    params = test_params(kSize=21, kDataFrameType="MQF")
    params.Hasher = kp.IntegerHasher(21)

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

        kf = kp.kDataFrameMQF(kSize)
        self.assertTrue(kf.empty())

        kf.insert(_kmer, rand_count)
        self.assertFalse(kf.empty())
        self.assertEqual(kf.count(_kmer), rand_count)

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
            c = kFrame.count(kmer[0])
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
                self.assertEqual(kFrames[i].count(kmer[0]), 0)

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
                count = it.getKmerCount()
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
                count = it.getKmerCount()
                kframe_kmers_counts.add(count)
                self.assertTrue(it.getHashedKmer() in kmers_hash_values)
                self.assertTrue(it.getKmerCount() in inserted_counts)
                it.next()

        shutil.rmtree(my_tmpdir)

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
