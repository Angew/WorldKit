import math
import unittest

import worldkit.placer as placer


class TestPoint(unittest.TestCase):
	def _testClose(self, isClose):
		# Exact/wide off
		self.assertTrue(isClose(placer.Point(3, 2), placer.Point(3, 2)))
		self.assertTrue(isClose(placer.Point(-1.0, 3.5), placer.Point(-1, 3.5)))
		self.assertFalse(isClose(placer.Point(-4, 0), placer.Point(0, -4)))
		self.assertFalse(isClose(placer.Point(2.2, 4.5), placer.Point(2.1, 4.5)))
		self.assertTrue(isClose(placer.Point(), placer.Point(0, 0)))
		# Close/almost close
		self.assertTrue(isClose(placer.Point(2, 1e-7), placer.Point(2, 2e-7)))
		self.assertTrue(isClose(placer.Point(-3, 2e-7), placer.Point(-3 + 1e-7, 3e-7)))
		self.assertFalse(isClose(placer.Point(3e-6, -2.6), placer.Point(5e-6, -2.6)))
		self.assertFalse(isClose(placer.Point(15.2, 6.3), placer.Point(15.2, 6.3 - 2e-6)))
		self.assertFalse(isClose(placer.Point(-14.4, -1e-6), placer.Point(-14.4 - 1.5e-6, 0.5e-6)))

	def test_close_points(self):
		self._testClose(lambda l, r: placer.Point.close(l, r))

	def test_close_scalars(self):
		self.assertTrue(placer.Point.close(3, 3))
		self.assertTrue(placer.Point.close(-4.7, -4.7))
		self.assertTrue(placer.Point.close(-4.7, -4.70000005))
		self.assertFalse(placer.Point.close(15.2, 15.200005))
		self.assertFalse(placer.Point.close(1, 2))

	def test_isClose(self):
		self._testClose(lambda l, r: l.isClose(r))

	def test_normalisation(self):
		norm_0_m3 = placer.Point(0, -3).normalisation()
		self.assertEqual(norm_0_m3[0], placer.Point(0, -1))
		self.assertEqual(norm_0_m3[1], 3)
		norm_1p5_2 = placer.Point(1.5, 2).normalisation()
		self.assertTrue(norm_1p5_2[0].isClose(placer.Point(0.6, 0.8)))
		self.assertAlmostEqual(norm_1p5_2[1], 2.5)
		norm_m1_0 = placer.Point(-1, 0).normalisation()
		self.assertEqual(norm_m1_0[0], placer.Point(-1, 0))
		self.assertEqual(norm_m1_0[1], 1)
		with self.assertRaises(ValueError):
			placer.Point().normalisation()
		with self.assertRaises(ValueError):
			placer.Point(0, 0).normalisation()
		norm_0_0 = placer.Point(0, 0).normalisation(onZero = placer.Point(42, 42))
		self.assertEqual(norm_0_0[0], placer.Point(42, 42))
		self.assertEqual(norm_0_0[1], 0)
		norm_1em8_2em8 = placer.Point(1e-8, 2e-8).normalisation(onZero = placer.Point(-6, 4))
		self.assertEqual(norm_1em8_2em8[0], placer.Point(-6, 4))
		self.assertEqual(norm_1em8_2em8[1], 0)

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

	def test__addition(self):
		self.assertEqual(placer.Point(1, 1) + placer.Point(0, 2), placer.Point(1, 3))
		self.assertEqual(placer.Point(-2.5, 2.5) + placer.Point(2.5, -2.5), placer.Point())
		
	def test__division(self):
		self.assertEqual(placer.Point(1, 4) / 4, placer.Point(0.25, 1))
		self.assertEqual(placer.Point(-3, 1.4) / 1, placer.Point(-3, 1.4))
		self.assertEqual(placer.Point(3.5, -2.2) / 0.25, placer.Point(14, -8.8))

	def test__equality(self):
		self.assertEqual(placer.Point(3, 2), placer.Point(3, 2))
		self.assertEqual(placer.Point(-1.0, 3.5), placer.Point(-1, 3.5))
		self.assertNotEqual(placer.Point(-4, 0), placer.Point(0, -4))
		self.assertNotEqual(placer.Point(2.2, 4.5), placer.Point(2.1, 4.5))
		self.assertEqual(placer.Point(), placer.Point(0, 0))

	def test__multiplication(self):
		self.assertEqual(placer.Point(2, 0) * 3.5, placer.Point(7, 0))
		self.assertEqual(placer.Point(4, -2.4) * 0, placer.Point())
		self.assertEqual(-1.7 * placer.Point(0, 2), placer.Point(0, -3.4))
		self.assertEqual(0 * placer.Point(-6.1, 14.42), placer.Point())
		
	def test__subtraction(self):
		self.assertEqual(placer.Point(4, 2) - placer.Point(1, 1), placer.Point(3, 1))
		self.assertEqual(placer.Point(-3.5, 13.2) - placer.Point(2.8, 0), placer.Point(-6.3, 13.2))
		self.assertEqual(placer.Point(1.2, -3.1) - placer.Point(1.2, -3.1), placer.Point())


if __name__ == '__main__':
	unittest.main()
