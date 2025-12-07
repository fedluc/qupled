from pathlib import Path

import sqlalchemy as sql

from .scheme_tables import SchemeTables
from .base_tables import ConflictMode, RunStatus

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
        engine = sql.create_engine(f"sqlite:///{database_path}")
        # Set sqlite properties
        DataBaseHandler._set_sqlite_pragma(engine)
        # Create tables
        self.table_metadata = sql.MetaData()
        self.scheme_tables = SchemeTables(engine)

    def insert_scheme_run(self, inputs):
        """ """
        self.scheme_tables.insert_run(inputs)

    def insert_scheme_inputs(self, inputs: dict[str, any]):
        """ """
        self.scheme_tables.insert_inputs(inputs)

    def insert_scheme_results(
        self,
        results: dict[str, any],
        conflict_mode: ConflictMode = ConflictMode.FAIL,
    ):
        """ """
        self.scheme_tables.insert_results(results, conflict_mode)

    def inspect_scheme_runs(self) -> list[dict[str, any]]:
        """ """
        return self.scheme_tables.inspect_runs()

    def update_scheme_run_status(self, status: RunStatus) -> None:
        """ """
        self.scheme_tables.update_run_status(status)

    def get_scheme_run(
        self,
        run_id: int,
        input_names: list[str] | None = None,
        result_names: list[str] | None = None,
    ) -> dict:
        """ """
        return self.scheme_tables.get_run(run_id, input_names, result_names)

    def get_scheme_inputs(self, run_id: int, names: list[str] | None = None) -> dict:
        """ """
        return self.scheme_tables.get_inputs(run_id, names)

    def get_scheme_results(self, run_id: int, names: list[str] | None = None) -> dict:
        """ """
        return self.scheme_tables.get_results(run_id, names)

    def delete_scheme_run(self, run_id: int) -> None:
        """ """
        self.scheme_tables.delete_run(run_id)

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
