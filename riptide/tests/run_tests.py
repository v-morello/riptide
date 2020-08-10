import os
import pytest


def test():
    """
    Run the test suite.
    """
    tests_dir = os.path.dirname(__file__)
    args = ['--verbose', tests_dir]
    pytest.main(args)