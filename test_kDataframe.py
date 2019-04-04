from __future__ import print_function
import kProcessor as kp
import unittest
import random

def generate_kmers(kSize, kmers_no):
    kmers_list = []
    for i in range(kmers_no):
        _kmer = "".join([random.choice('ACGT') for _ in range(kSize)])
        _count = random.randint(1, 1000)
        kmers_list.append((_kmer, _count))

    return kmers_list

class TestKDataFrame(unittest.TestCase):

    fastqFiles = ["testData/test.noN.fastq", "testData/test2.noN.fastq", "testData/test2.noN.fastq"]

    def test_emptykDataFrame(self):
        print(self._testMethodName)

        framesToBeTested = []
        kSizes = [21, 31]
        for kSize in kSizes:
            framesToBeTested.append(kp.kDataFrameMQF(kSize))

        for kSize in kSizes:
            framesToBeTested.append(kp.kDataFrameMAP(kSize))

        for kFrame in framesToBeTested:
            self.assertTrue(kFrame.empty())

    def test_build_kDataFrames(self):
        k = 31
        framesToBeTested = [[],[]]

        for file in self.fastqFiles:
            framesToBeTested[0].append(kp.kDataFrameMQF(31))
            kp.parseSequences(file, 1, framesToBeTested[0][-1])
            framesToBeTested[1].append(kp.kDataFrameMAP(31))
            kp.parseSequences(file, 1, framesToBeTested[1][-1])

        self.framesToBeTested = framesToBeTested

    def test_insertOneTime(self):
        print(self._testMethodName)
        kSize = 31
        _kmer = "".join([random.choice('ACGT') for _ in range(kSize)])
        rand_count = random.randint(1, 1000)

        df_mqf = kp.kDataFrameMQF(kSize)
        self.assertTrue(df_mqf.empty())

        df_mqf.insert(_kmer, rand_count)
        self.assertFalse(df_mqf.empty())
        self.assertEqual(df_mqf.count(_kmer), rand_count)

        df_map = kp.kDataFrameMAP(kSize)
        self.assertTrue(df_map.empty())
        
        df_map.insert(_kmer, rand_count)
        self.assertFalse(df_map.empty())
        self.assertEqual(df_map.count(_kmer), rand_count)

    def test_insertNTimes(self):
        print(self._testMethodName)
        kSize = 31

        # Create random kmers list
        kmers_list = generate_kmers(kSize, kmers_no=20)

        # Empty kDataFrames
        kFrames = [kp.kDataFrameMQF(kSize), kp.kDataFrameMAP(kSize)]
        self.assertTrue(kFrames[0].empty())
        self.assertTrue(kFrames[1].empty())

        # Count inserted kmers for [kDataFrameMQF, kDataFrameMAP]
        insertedKmers = [0, 0]

        # Insert all kmers
        for i in range(len(kFrames)):
            for kmer in kmers_list:
                if(kFrames[i].insert(kmer[0], kmer[1])):
                    insertedKmers[i] += 1

        # Assert all kmers are inserted
        self.assertEqual(insertedKmers[0], len(kmers_list))  # MQF
        self.assertEqual(insertedKmers[1], len(kmers_list))  # MAP

        # Verify all inserted kmers
        for i in range(len(kFrames)):
            for kmer in kmers_list:
                c = kFrames[i].count(kmer[0])
                self.assertEqual(c, kmer[1])

    def test_eraseKmers(self):
        print(self._testMethodName)
        kSize = 31

        # Create random kmers list
        kmers_list = generate_kmers(kSize, kmers_no=20)

        # Empty kDataFrames
        kFrames = [kp.kDataFrameMQF(kSize), kp.kDataFrameMAP(kSize)]
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
        kSize = 31

        # Create random kmers list
        kmers_list = generate_kmers(kSize, kmers_no=20)

        # Empty kDataFrames
        kFrames = [kp.kDataFrameMQF(kSize), kp.kDataFrameMAP(kSize)]
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
            while(it != kFrame.end()):
                count = it.getKmerCount()
                kframe_kmers_counts.add(count)
                it.next()

            self.assertEqual(len(kframe_kmers_counts.intersection(inserted_counts)), len(inserted_counts))

    def test_saveAndIterateOverAllKmers(self):
        print(self._testMethodName)
        kSize = 31

        # Create random kmers list
        kmers_list = generate_kmers(kSize, kmers_no=20)

        # Empty kDataFrames
        kFrames = [kp.kDataFrameMQF(kSize), kp.kDataFrameMAP(kSize)]
        self.assertTrue(kFrames[0].empty())
        self.assertTrue(kFrames[1].empty())

        # Insert kmers
        for i in range(len(kFrames)):
            for kmer in kmers_list:
                self.assertTrue(kFrames[i].insert(kmer[0], kmer[1]))

        # Get all inserted kmers counts
        inserted_counts = set([kmer[1] for kmer in kmers_list])

        fileNames = ["tmp.kdataframe.mqf", "tmp.kdataframe.map"]

        for i in range(len(kFrames)):
            kFrames[i].save(fileNames[i])

        loaded_kFrames = []
        for i in range(len(kFrames)):
            loaded_kFrames.append(kp.kDataFrame.load(fileNames[i]))

        for kFrame in loaded_kFrames:
            it = kFrame.begin()
            kframe_kmers_counts = set()
            while(it != kFrame.end()):
                count = it.getKmerCount()
                kframe_kmers_counts.add(count)
                it.next()

            self.assertEqual(len(kframe_kmers_counts.intersection(inserted_counts)), len(inserted_counts))

    def test_unioinTest(self):
        print(self._testMethodName)
        kSize = 31

        # Empty kDataFrames
        ## auto binding to vector<kDataFrame*>
        kFrames_vec = [kp.kDataFrameMQF(kSize), kp.kDataFrameMAP(kSize)]
        self.assertTrue(kFrames_vec[0].empty())
        self.assertTrue(kFrames_vec[1].empty())

        # Create random kmers list
        kmers_list1 = generate_kmers(kSize, kmers_no=20)
        inserted_counts1 = set([kmer[1] for kmer in kmers_list1])

        kmers_list2 = generate_kmers(kSize, kmers_no=20)
        inserted_counts2 = set([kmer[1] for kmer in kmers_list2])

        # Inserting Kmers

        for kmer in kmers_list1:
            self.assertTrue(kFrames_vec[0].insert(kmer[0], kmer[1]))

        for kmer in kmers_list2:
            self.assertTrue(kFrames_vec[1].insert(kmer[0], kmer[1]))

        # Apply kFrameUnion
        union_kFrame = kp.kFrameUnion(kFrames_vec)

        # Total kmers extracted from union kFrames_vec
        kframe_kmers_counts = set()

        it = union_kFrame.begin()

        while(it != union_kFrame.end()):
            count = it.getKmerCount()
            kframe_kmers_counts.add(count)
            it.next()

        self.assertEqual(len(kframe_kmers_counts), len(inserted_counts1.union(inserted_counts2)))

    def test_intersectTest(self):
        print(self._testMethodName)
        kSize = 31

        # Empty kDataFrames
        ## auto binding to vector<kDataFrame*>
        kFrames_vec = [kp.kDataFrameMQF(kSize), kp.kDataFrameMQF(kSize)]

        self.assertTrue(kFrames_vec[0].empty())
        self.assertTrue(kFrames_vec[1].empty())

        # Create random kmers list
        kmers_list1 = generate_kmers(kSize, kmers_no=20)
        kmers_list2 = generate_kmers(kSize, kmers_no=10)

        # Replicate some kmers from kmers_list1 in kmers_list2 to make sure len(intersection) > 0
        kmers_list2 += kmers_list1[0:10]

        # Inserting Kmers
        for kmer in kmers_list1:
            self.assertTrue(kFrames_vec[0].insert(kmer[0], kmer[1]))

        for kmer in kmers_list2:
            self.assertTrue(kFrames_vec[1].insert(kmer[0], kmer[1]))

        # Apply kFrameIntersect
        intersect_kFrame = kp.kFrameIntersect(kFrames_vec)

        self.assertGreater(intersect_kFrame.size(), 1)


    def test_differenceTest(self):
        print(self._testMethodName)
        kSize = 31

        # Empty kDataFrames
        ## auto binding to vector<kDataFrame*>
        kFrames_vec = [kp.kDataFrameMQF(kSize), kp.kDataFrameMQF(kSize)]

        self.assertTrue(kFrames_vec[0].empty())
        self.assertTrue(kFrames_vec[1].empty())

        # Create random kmers list
        kmers_list1 = generate_kmers(kSize, kmers_no=20)
        kmers_list2 = generate_kmers(kSize, kmers_no=10)

        # Replicate some kmers from kmers_list1 in kmers_list2 to make sure len(intersection) > 0
        kmers_list2 += kmers_list1[0:10]

        # Inserting Kmers
        for kmer in kmers_list1:
            self.assertTrue(kFrames_vec[0].insert(kmer[0], kmer[1]))

        for kmer in kmers_list2:
            self.assertTrue(kFrames_vec[1].insert(kmer[0], kmer[1]))

        # Apply kFrameIntersect
        diff_kFrame = kp.kFrameIntersect(kFrames_vec)

        self.assertGreater(diff_kFrame.size(), 1)



if __name__ == '__main__':
    unittest.main(verbosity=3)
