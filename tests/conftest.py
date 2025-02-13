from threading import Thread

import pytest
from dotenv import load_dotenv

from app import create_app

load_dotenv()
# Load secrets from secrets.env file if available
load_dotenv("secrets.env")


@pytest.fixture(scope='session')
def app():
    # Create app with shutdown hook
    app = create_app()

    # Hack: Start Flask Werkzeug server as a daemon thread and let it dangle (skip `thread.join()`). When the tests
    # finish, it will be the only thread left running and the Python process exits when only daemon threads are left.
    # Details here: https://docs.python.org/3/library/threading.html#thread-objects
    thread = Thread(target=app.run, daemon=True, kwargs=dict(port=5000))
    thread.start()

    yield app


@pytest.fixture(scope='module')
def client(app):
    return app.test_client()
