import pytest
from tests.utilities import *
from app.common import *


"""
Unit Tests
----------
"""

@pytest.mark.parametrize("input, expected_result", [
    ("", False),
    ("1", True),
    ("1.2", False),
    ("432", True),
    ("True", False)
])
def test_is_int(input, expected_result):
    assert is_int(input) == expected_result
