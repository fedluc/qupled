import os
import pytest
import sqlalchemy as sql
from sqlalchemy import inspect
from qupled.database import DataBaseHandler
from sqlalchemy import event
from sqlalchemy.engine import Engine
from sqlalchemy.exc import OperationalError


def test_database_handler_initialization_with_default_name():
    handler = DataBaseHandler()
    assert handler.database_name == DataBaseHandler.DEFAULT_DATABASE_NAME
    assert handler.engine.url.database == DataBaseHandler.DEFAULT_DATABASE_NAME
    assert handler.run_id is None
    inspector = inspect(handler.engine)
    assert set(inspector.get_table_names()) == {
        DataBaseHandler.RUNS_TABLE_NAME,
        DataBaseHandler.INPUTS_TABLE_NAME,
        DataBaseHandler.RESULTS_TABLE_NAME,
    }


def test_database_handler_initialization_with_custom_name():
    database_name = "custom.db"
    handler = DataBaseHandler(database_name="custom.db")
    assert handler.database_name == database_name
    assert handler.engine.url.database == database_name
    assert handler.run_id is None
    inspector = inspect(handler.engine)
    assert set(inspector.get_table_names()) == {
        DataBaseHandler.RUNS_TABLE_NAME,
        DataBaseHandler.INPUTS_TABLE_NAME,
        DataBaseHandler.RESULTS_TABLE_NAME,
    }


def test_set_sqlite_pragma_enforces_foreign_keys():
    engine = sql.create_engine("sqlite:///:memory:")
    DataBaseHandler._set_sqlite_pragma(engine)
    with engine.connect() as connection:
        result = connection.execute(sql.text("PRAGMA foreign_keys")).fetchone()
        assert result[0] == 1


def test_set_sqlite_pragma_does_not_raise_error_on_invalid_engine():
    invalid_engine = None
    with pytest.raises(sql.exc.InvalidRequestError):
        DataBaseHandler._set_sqlite_pragma(invalid_engine)
