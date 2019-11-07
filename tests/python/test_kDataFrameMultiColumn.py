import os
import shutil
import unittest
import random
import kProcessor as kp
import random
import sys
from params import test_params

class TestkDataFrameMultiColumn(unittest.TestCase):
    
    def main_test(self, kf_type):
        params = test_params(kSize=21, kDataFrameType=kf_type)

        kmers_list = params.generate_kmers(kmers_no=1000)

        kFrame = params.create_new_kf()

        # Insert all kmers
        for kmer in kmers_list:
            kFrame.insert(kmer[0], kmer[1])

        
        kFrame.addColumn_int("intColumn")
        kFrame.addColumn_bool("boolColumn")
        kFrame.addColumn_double("doubleColumn")

        simColumns = dict()
        it = kFrame.begin()

        while it != kFrame.end():
            kmer = it.getKmer()
            _double = random.uniform(1.0,100.0)
            _bool = bool(random.getrandbits(1))
            _int = random.randint(1,1000)
            simColumns[kmer] = (_int, _bool, _double)
            kFrame.setKmerColumnValue_int("intColumn",kmer,_int)
            kFrame.setKmerColumnValue_bool("boolColumn",kmer,_bool)
            kFrame.setKmerColumnValue_double("doubleColumn",kmer,_double)
            it.next()


        for kmer, values in simColumns.items():
            _int = kFrame.getKmerColumnValue_int("intColumn",kmer)
            _bool = kFrame.getKmerColumnValue_bool("boolColumn",kmer)
            _double = kFrame.getKmerColumnValue_double("doubleColumn",kmer)
            

            golden_int = values[0]
            golden_bool = values[1]
            golden_double = values[2]

            self.assertEqual(golden_int,_int)
            self.assertEqual(golden_bool,_bool)
            self.assertEqual(golden_double,_double)


    def test_MultiKF_MQF(self):
        self.main_test("MQF")

    # def test_MultiKF_MAP(self):
    #     self.main_test("MAP")

    # def test_MultiKF_PHMAP(self):
    #     self.main_test("PHMAP")