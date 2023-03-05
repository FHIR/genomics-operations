import os

import pytest

from app import create_app


@pytest.fixture(scope='module')
def client():
    with create_app().test_client() as c:
        yield c
