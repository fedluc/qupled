import os
from pathlib import Path

import pytest
import sqlalchemy as sql

from qupled.database.database_handler import (
    DataBaseHandler,
    FiniteSizeCorrectionTables,
    SchemeTables,
    BLOB_STORAGE_DIRECTORY,
    DATABASE_DIRECTORY,
    DEFAULT_DATABASE_NAME,
)


@pytest.fixture
def db_handler():
    handler = DataBaseHandler()
    yield handler
    database_url = handler.engine.url.database
    if os.path.exists(database_url):
        os.remove(database_url)


def test_database_handler_initialization_with_default_name(db_handler):
    expected_database_path = Path(DATABASE_DIRECTORY) / DEFAULT_DATABASE_NAME
    expected_blob_storage = (
        Path(DATABASE_DIRECTORY) / BLOB_STORAGE_DIRECTORY / DEFAULT_DATABASE_NAME
    )
    assert db_handler.blob_storage == str(expected_blob_storage)
    assert db_handler.engine.url.database == str(expected_database_path)
    assert isinstance(db_handler.scheme_tables, SchemeTables)
    assert isinstance(db_handler.fsc_tables, FiniteSizeCorrectionTables)


def test_database_handler_initialization_with_custom_name():
    database_name = "custom.db"
    expected_database_path = Path(DATABASE_DIRECTORY) / database_name
    expected_blob_storage = (
        Path(DATABASE_DIRECTORY) / BLOB_STORAGE_DIRECTORY / database_name
    )
    db_handler = DataBaseHandler(database_name="custom.db")
    assert db_handler.blob_storage == str(expected_blob_storage)
    assert db_handler.engine.url.database == str(expected_database_path)
    assert isinstance(db_handler.scheme_tables, SchemeTables)
    assert isinstance(db_handler.fsc_tables, FiniteSizeCorrectionTables)
    if os.path.exists(expected_database_path):
        os.remove(expected_database_path)


def test_set_sqlite_pragma_valid_engine():
    engine = sql.create_engine("sqlite:///:memory:")
    DataBaseHandler._set_sqlite_pragma(engine)
    with engine.connect() as connection:
        result = connection.execute(sql.text("PRAGMA foreign_keys")).fetchone()
        assert result[0] == 1


def test_set_sqlite_pragma_invalid_engine():
    invalid_engine = None
    with pytest.raises(sql.exc.InvalidRequestError):
        DataBaseHandler._set_sqlite_pragma(invalid_engine)
