import unittest
import aipy.dsp as dsp

class TestWindows(unittest.TestCase):
    def test_kaiser2(self):
        win = dsp.gen_window(1024, 'kaiser2')
        self.assertAlmostEqual(win[0], 0.01147993)
        self.assertAlmostEqual(win[1], 0.01192681)
        self.assertAlmostEqual(win[2], 0.01238142)

if __name__ == '__main__':
    unittest.main()
