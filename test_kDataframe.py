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

    def test_instances_types(self):
        """
        check the types of the variables
        """
        print(self._testMethodName)

        _int_vector = kp.IntVector()
        self.assertIsInstance(_int_vector, kp.IntVector)

        _kDataFrameMQF = kp.kDataFrameMQF(31)
        self.assertIsInstance(_kDataFrameMQF, kp.kDataFrameMQF)

        _kDataFrameMAP = kp.kDataFrameMAP(31)
        self.assertIsInstance(_kDataFrameMAP, kp.kDataFrameMAP)

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

        # Count inserted kmers for [kDataFrameMQF, kDataFrameMAP]
        insertedKmers = [0, 0]

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

        verifiedKmers = [0, 0]
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

    def test_insertNTimes(self):
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

    



if __name__ == '__main__':
    unittest.main(verbosity=3)
