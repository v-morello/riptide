import os
import pytest


def test():
    """
    Run the test suite.
    """
    tests_dir = os.path.join(os.path.dirname(__file__), '..', 'tests')
    args = ['--verbose', tests_dir]
    pytest.main(args)