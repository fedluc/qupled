from pathlib import Path

import sqlalchemy as sql

from qupled.database.scheme_tables import SchemeTables
from qupled.database.finite_size_correction_tables import FiniteSizeCorrectionTables

BLOB_STORAGE_DIRECTORY = "blob_data"
DATABASE_DIRECTORY = "qupled_store"
DEFAULT_DATABASE_NAME = "qupled.db"


class DataBaseHandler:
    """
    DataBaseHandler is a class for managing a SQLite database that stores information
    about runs, inputs, and results. It provides methods for inserting, retrieving,
    and deleting data, as well as managing the database schema."
    """

    def __init__(self, database_name: str | None = None):
        """ """
        # Database path
        database_name = (
            DEFAULT_DATABASE_NAME if database_name is None else database_name
        )
        database_path = Path(DATABASE_DIRECTORY) / database_name
        database_path.parent.mkdir(parents=True, exist_ok=True)
        # Blob data storage
        self.blob_storage = (
            Path(DATABASE_DIRECTORY) / BLOB_STORAGE_DIRECTORY / database_name
        )
        self.blob_storage.mkdir(parents=True, exist_ok=True)
        self.blob_storage = str(self.blob_storage)
        # Create database
        self.engine = sql.create_engine(f"sqlite:///{database_path}")
        # Set sqlite properties
        DataBaseHandler._set_sqlite_pragma(self.engine)
        # Create tables
        self.scheme_tables = SchemeTables(self.engine)
        self.fsc_tables = FiniteSizeCorrectionTables(self.engine)

    @staticmethod
    def _set_sqlite_pragma(engine):
        """
        Configures the SQLite database engine to enforce foreign key constraints.

        This function sets up a listener for the "connect" event on the provided
        SQLAlchemy engine. When a new database connection is established, it executes
        the SQLite PRAGMA statement to enable foreign key support.

        Args:
            engine (sqlalchemy.engine.Engine): The SQLAlchemy engine instance to configure.

        Notes:
            SQLite does not enforce foreign key constraints by default. This function
            ensures that foreign key constraints are enabled for all connections made
            through the provided engine.
        """

        @sql.event.listens_for(engine, "connect")
        def _set_pragma(dbapi_connection, connection_record):
            cursor = dbapi_connection.cursor()
            cursor.execute("PRAGMA foreign_keys=ON")
            cursor.execute("PRAGMA journal_mode=WAL")
            cursor.close()
