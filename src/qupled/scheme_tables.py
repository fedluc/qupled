from datetime import datetime
import sqlalchemy as sql

from . import native
from . import base_tables as base
from . import hf

INPUT_TABLE_NAME = "inputs"
RESULT_TABLE_NAME = "results"
RUN_TABLE_NAME = "runs"


class SchemeTables(base.BaseTables):
    """ """

    def __init__(self, engine: sql.Engine):
        """
        Initializes the DataBaseHandler instance.

        Args:
            database_name (str | None, optional): The name of the database file. If not provided,
                the default database name (`DEFAULT_DATABASE_NAME`) will be used.

        Attributes:
            database_name (str): The name of the database file being used.
            engine (sqlalchemy.engine.Engine): The SQLAlchemy engine connected to the SQLite database.
            table_metadata (sqlalchemy.MetaData): Metadata object for managing table schemas.
            run_table (sqlalchemy.Table): The table schema for storing run information.
            input_table (sqlalchemy.Table): The table schema for storing input data.
            result_table (sqlalchemy.Table): The table schema for storing result data.
            run_id (int | None): The ID of the current run, or None if no run is active.
        """
        super().__init__(
            self, engine, RUN_TABLE_NAME, INPUT_TABLE_NAME, RESULT_TABLE_NAME
        )
        self._build_tables()

    def delete_run(self, run_id: int) -> None:
        """
        Deletes a run entry from the database based on the provided run ID.

        Args:
            run_id (int): The unique identifier of the run to be deleted.

        Returns:
            None
        """
        native.delete_blob_data_on_disk(self.engine.url.database, run_id)
        super().delete_run(run_id)

    def _build_run_table(self):
        """
        Builds the SQLAlchemy table object for the "runs" table in the database.

        This method defines the schema for the "runs" table, including its columns,
        data types, constraints, and metadata. The table includes the following columns:

        - PRIMARY_KEY: An auto-incrementing integer that serves as the primary key.
        - THEORY: A string representing the theory associated with the run (non-nullable).
        - COUPLING: A float representing the coupling value (non-nullable).
        - DEGENERACY: A float representing the degeneracy value (non-nullable).
        - DATE: A string representing the date of the run (non-nullable).
        - TIME: A string representing the time of the run (non-nullable).
        - STATUS: A string representing the status of the run (non-nullable).

        After defining the table schema, the method creates the table in the database
        using the `_create_table` method.

        Returns:
            sqlalchemy.Table: The constructed SQLAlchemy table object for the "runs" table.
        """
        table = sql.Table(
            self.run_table_name,
            self.table_metadata,
            sql.Column(
                base.TableKeys.PRIMARY_KEY.value,
                sql.Integer,
                primary_key=True,
                autoincrement=True,
            ),
            sql.Column(
                base.TableKeys.THEORY.value,
                sql.String,
                nullable=False,
            ),
            sql.Column(
                base.TableKeys.COUPLING.value,
                sql.Float,
                nullable=False,
            ),
            sql.Column(
                base.TableKeys.DEGENERACY.value,
                sql.Float,
                nullable=False,
            ),
            sql.Column(
                base.TableKeys.DATE.value,
                sql.String,
                nullable=False,
            ),
            sql.Column(
                base.TableKeys.TIME.value,
                sql.String,
                nullable=False,
            ),
            sql.Column(base.TableKeys.STATUS.value, sql.String, nullable=False),
        )
        self._create_table(table)
        return table

    def _insert_run(self, inputs: hf.Input, status: base.RunStatus):
        """
        Inserts a new run entry into the database.

        Args:
            inputs (any): An object containing the input data for the run.
                          Expected attributes include:
                          - theory: Theoretical data to be serialized into JSON.
                          - coupling: Coupling data to be serialized into JSON.
                          - degeneracy: Degeneracy data to be serialized into JSON.

        Side Effects:
            - Updates the `self.run_id` attribute with the primary key of the newly inserted run.

        Notes:
            - The current date and time are automatically added to the entry.
            - The input data is serialized into JSON format before insertion.
        """
        now = datetime.now()
        data = {
            base.TableKeys.THEORY.value: inputs.theory,
            base.TableKeys.COUPLING.value: inputs.coupling,
            base.TableKeys.DEGENERACY.value: inputs.degeneracy,
            base.TableKeys.DATE.value: now.date().isoformat(),
            base.TableKeys.TIME.value: now.time().isoformat(),
            base.TableKeys.STATUS.value: status.value,
        }
        statement = sql.insert(self.run_table).values(data)
        result = self._execute(statement)
        if run_id := result.inserted_primary_key:
            self.run_id = run_id[0]
