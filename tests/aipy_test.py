# -*- coding: utf-8 -*-

"""Unit test suite for the AIPY package."""

import unittest

import _alm_test
import _healpix_test
import amp_test
import coord_test
import deconv_test
import miriad_test
import phs_test
import phs_benchmark
import scripting_test

class TestSuite(unittest.TestSuite):
        """A unittest.TestSuite class which contains all of the package unit tests."""

        def __init__(self):
                """Setup the AIPY package unit test suite."""

                unittest.TestSuite.__init__(self)

                self.addTest(_alm_test.TestSuite())
                self.addTest(_healpix_test.TestSuite())
                self.addTest(amp_test.TestSuite())
                self.addTest(coord_test.TestSuite())
                self.addTest(deconv_test.TestSuite())
                self.addTest(miriad_test.TestSuite())
                self.addTest(phs_test.TestSuite())
                self.addTest(phs_benchmark.TestSuite())
                self.addTest(scripting_test.TestSuite())

def main(opts=None, args=None):
    """Function to call all of the lsl tests."""

    if opts is not None:
        if opts.verbose:
            level = 2
        else:
            level = 1
    else:
        level = 2

    suite = TestSuite()
    runner = unittest.TextTestRunner(verbosity = level)
    runner.run(suite)


if __name__  == '__main__':
    import optparse

    parser = optparse.OptionParser(usage = "python %prog [options]", description = __doc__)
    parser.add_option("-v", "--verbose", action = "store_true", dest = "verbose", default = False,
                      help = "extra print output")
    (opts, args) = parser.parse_args()

    main(opts, args)
