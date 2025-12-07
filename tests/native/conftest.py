import os
import pytest
from qupled.database import DEFAULT_DATABASE_NAME


@pytest.fixture(autouse=True)
def run_after_each_test():
    yield
    database_name = DEFAULT_DATABASE_NAME
    if os.path.exists(database_name):
        os.remove(database_name)
