"""
Run all test cases.
"""

import os
import unittest

import reflexible as rf


def suite():
    this_dir = os.path.dirname(__file__)
    return unittest.TestLoader().discover(
        start_dir=this_dir, pattern="test_*.py")


def test():
    """
    test()

    Run all the tests in the test suite.
    """
    rf.print_versions()
    return unittest.TextTestRunner().run(suite())


if __name__ == '__main__':
    test()
