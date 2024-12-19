import pytest
from dotenv import load_dotenv

from app import create_app

load_dotenv()
# Load secrets from secrets.env file if available
load_dotenv("secrets.env")


@pytest.fixture(scope='module')
def client():
    with create_app().test_client() as c:
        yield c
