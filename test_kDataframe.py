import kProcessor as kp
import unittest

class TestKDataFrame(unittest.TestCase):

    def test_instances_types(self):
        print(self._testMethodName)

        _int_vector = kp.IntVector()
        self.assertIsInstance(_int_vector, kp.IntVector)

        _kDataFrameMQF = kp.kDataFrameMQF(31)
        self.assertIsInstance(_kDataFrameMQF, kp.kDataFrameMQF)

        _kDataFrameMAP = kp.kDataFrameMAP(31)
        self.assertIsInstance(_kDataFrameMAP, kp.kDataFrameMAP)


if __name__ == '__main__':
    unittest.main()