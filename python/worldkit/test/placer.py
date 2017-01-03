import math
import unittest

import worldkit.placer as placer


class TestPoint(unittest.TestCase):
	def test_size(self):
		self.assertAlmostEqual(placer.Point(-2, 0).size(), 2.0)
		self.assertAlmostEqual(placer.Point(1, 1).size(), math.sqrt(2))
		self.assertEqual(placer.Point().size(), 0)
		self.assertEqual(placer.Point(0, 0).size(), 0)

	def test_size2(self):
		self.assertEqual(placer.Point(-2, 0).size2(), 4.0)
		self.assertEqual(placer.Point(3, 2).size2(), 13.0)
		self.assertEqual(placer.Point().size2(), 0)
		self.assertEqual(placer.Point(0, 0).size2(), 0)

	def test_normalisation(self):
		norm_3_0 = placer.Point(0, -3).normalisation()
		self.assertEqual(norm_3_0[0], placer.Point(0, -1))


if __name__ == '__main__':
	unittest.main()
