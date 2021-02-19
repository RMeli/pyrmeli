"""
Unit and regression test for the pyrmeli package.
"""

# Import package, test suite, and other packages as needed
import pyrmeli
import pytest
import sys

def test_pyrmeli_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "pyrmeli" in sys.modules
