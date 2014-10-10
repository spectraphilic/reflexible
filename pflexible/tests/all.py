"""
Run all test cases.
"""

import sys
import os
import unittest

import pflexible as pf


def suite():
    this_dir = os.path.dirname(__file__)
    return unittest.TestLoader().discover(
        start_dir=this_dir, pattern="test_*.py")


def test(only_versions=False):
    """
    test()

    Run all the tests in the test suite.
    """
    pf.print_versions()
    if not only_versions:
        return unittest.TextTestRunner().run(suite())


if __name__ == '__main__':

    # Handle some global flags (i.e. only useful for test_all.py)
    only_versions = False
    args = sys.argv[:]
    for arg in args:
        if arg in ['-print-versions']:
            only_versions = True
            sys.argv.remove(arg)
    test(only_versions)
