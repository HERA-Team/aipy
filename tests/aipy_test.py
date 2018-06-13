# -*- coding: utf-8 -*-

"""
Unit test suite for the AIPY package.
"""

from __future__ import print_function, division, absolute_import

import unittest

from . import _alm_test
from . import _healpix_test
from . import amp_test
from . import coord_test
from . import deconv_test
from . import miriad_test
from . import optimize_test
from . import phs_test
from . import phs_benchmark
from . import scripting_test
from . import twodgauss_test

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
                self.addTest(optimize_test.TestSuite())
                self.addTest(phs_test.TestSuite())
                self.addTest(phs_benchmark.TestSuite())
                self.addTest(scripting_test.TestSuite())
                self.addTest(twodgauss_test.TestSuite())

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
