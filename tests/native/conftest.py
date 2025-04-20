import os
import pytest
from qupled.database import DataBaseHandler


@pytest.fixture(autouse=True)
def run_after_each_test():
    yield
    database_name = DataBaseHandler.DEFAULT_DATABASE_NAME
    if os.path.exists(database_name):
        os.remove(database_name)
