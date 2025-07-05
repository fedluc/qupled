import io
import json
import os
import struct
from datetime import datetime
from enum import Enum
from collections.abc import Callable
from pathlib import Path

import numpy as np
import sqlalchemy as sql
from sqlalchemy.exc import NoSuchTableError
from sqlalchemy.dialects.sqlite import insert as sqlite_insert
import blosc2


class DataBaseHandler:
    """
    DataBaseHandler is a class for managing a SQLite database that stores information
    about runs, inputs, and results. It provides methods for inserting, retrieving,
    and deleting data, as well as managing the database schema."
    """

    BLOB_STORAGE_DIRECTORY = "blob_data"
    DATABASE_DIRECTORY = "qupled_store"
    DEFAULT_DATABASE_NAME = "qupled.db"
    FIXED_TABLE_NAME = "fixed"
    INPUT_TABLE_NAME = "inputs"
    RESULT_TABLE_NAME = "results"
    RUN_TABLE_NAME = "runs"

    class TableKeys(Enum):
        COUPLING = "coupling"
        DATE = "date"
        DEGENERACY = "degeneracy"
        NAME = "name"
        PRIMARY_KEY = "id"
        RUN_ID = "run_id"
        STATUS = "status"
        THEORY = "theory"
        TIME = "time"
        VALUE = "value"

    class RunStatus(Enum):
        RUNNING = "STARTED"
        SUCCESS = "SUCCESS"
        FAILED = "FAILED"

    class ConflictMode(Enum):
        FAIL = "FAIL"
        UPDATE = "UPDATE"

    def __init__(self, database_name: str | None = None):
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
        # Database path
        database_name = (
            self.DEFAULT_DATABASE_NAME if database_name is None else database_name
        )
        database_path = Path(self.DATABASE_DIRECTORY) / database_name
        database_path.parent.mkdir(parents=True, exist_ok=True)
        # Blob data storage
        self.blob_storage = (
            Path(self.DATABASE_DIRECTORY) / self.BLOB_STORAGE_DIRECTORY / database_name
        )
        self.blob_storage.mkdir(parents=True, exist_ok=True)
        self.blob_storage = str(self.blob_storage)
        # Create database
        self.engine = sql.create_engine(f"sqlite:///{database_path}")
        # Set sqlite properties
        DataBaseHandler._set_sqlite_pragma(self.engine)
        # Create tables
        self.table_metadata = sql.MetaData()
        self.run_table = self._build_run_table()
        self.input_table = self._build_inputs_table()
        self.result_table = self._build_results_table()
        self.run_id: int | None = None

    def insert_run(self, inputs):
        """
        Inserts a new run into the database by storing the provided inputs and results.

        Args:
            inputs (object): An object containing the input data for the run.
                             The attributes of this object will be converted to a dictionary.
            results (object): An object containing the result data for the run.
                              The attributes of this object will be converted to a dictionary.

        """
        self._insert_run(inputs, self.RunStatus.RUNNING)
        self.insert_inputs(inputs.__dict__)

    def insert_inputs(self, inputs: dict[str, any]):
        """
        Inserts input data into the database for the current run.

        Args:
            inputs (dict[str, any]): A dictionary containing input data to be inserted.
                The keys represent column names, and the values represent the data
                to be stored in the corresponding columns.

        Raises:
            ValueError: If `run_id` is None, indicating that no run is currently active.

        Notes:
            - The input data is serialized to JSON format before being inserted into
              the database.
            - The insertion is performed using the `_insert_from_dict` method, which
              maps the input values using the `sql_mapping` function.
        """
        if self.run_id is not None:
            sql_mapping = lambda value: (self._to_json(value))
            self._insert_from_dict(self.input_table, inputs, sql_mapping)

    def insert_results(
        self,
        results: dict[str, any],
        conflict_mode: ConflictMode = ConflictMode.FAIL,
    ):
        """
        Inserts the given results into the database table associated with this instance.

        Args:
            results (dict[str, any]): A dictionary where the keys are column names
                and the values are the corresponding data to be inserted.

        Notes:
            - This method requires that `self.run_id` is not None; otherwise, no insertion will occur.
            - The values in the `results` dictionary are converted to bytes using the `_to_bytes` method
              before being inserted into the database.
        """
        if self.run_id is not None:
            sql_mapping = lambda value: (self._to_bytes(value))
            self._insert_from_dict(
                self.result_table, results, sql_mapping, conflict_mode
            )

    def inspect_runs(self) -> list[dict[str, any]]:
        """
        Retrieve and inspect all runs from the database.

        This method executes a SQL SELECT statement on the `run_table` and retrieves
        all rows. Each row is converted into a dictionary where the keys are the column
        names and the values are the corresponding data.

        Returns:
            list[dict[str, any]]: A list of dictionaries, each representing a row
            from the `run_table`. The keys in the dictionary correspond to the column
            names, and the values are the respective data for each column.
        """
        statement = sql.select(self.run_table)
        rows = self._execute(statement).mappings().all()
        return [{key: row[key] for key in row.keys()} for row in rows]

    def update_run_status(self, status: RunStatus) -> None:
        """
        Update the status of a run in the database.

        Args:
            status (RunStatus): The new status to set for the run.

        Returns:
            None

        Notes:
            This method updates the status of the run identified by `self.run_id` in the run table.
            If `self.run_id` is None, no update is performed.
        """
        if self.run_id is not None:
            statement = (
                sql.update(self.run_table)
                .where(
                    self.run_table.c[self.TableKeys.PRIMARY_KEY.value] == self.run_id
                )
                .values({self.TableKeys.STATUS.value: status.value})
            )
            self._execute(statement)

    def get_run(
        self,
        run_id: int,
        input_names: list[str] | None = None,
        result_names: list[str] | None = None,
    ) -> dict:
        """
        Retrieve a run's data, including its inputs and results, from the database.

        Args:
            run_id (int): The unique identifier of the run to retrieve.
            input_names (list[str] | None): A list of input names to filter the inputs
                associated with the run. If None, all inputs are retrieved.
            result_names (list[str] | None): A list of result names to filter the results
                associated with the run. If None, all results are retrieved.

        Returns:
            dict: A dictionary containing the run's data, inputs, and results. The structure is:
                {
                    "RUN_TABLE_NAME": {<run_data>},
                    "INPUT_TABLE_NAME": [<inputs>],
                    "RESULT_TABLE_NAME": [<results>]
                If the run is not found, an empty dictionary is returned.
        """
        statement = sql.select(self.run_table).where(
            self.run_table.c[self.TableKeys.PRIMARY_KEY.value] == run_id
        )
        result = self._execute(statement).mappings().first()
        if result is not None:
            run_data = {key: result[key] for key in result.keys()}
            inputs = self.get_inputs(run_id, names=input_names)
            results = self.get_results(run_id, names=result_names)
            return {
                self.RUN_TABLE_NAME: run_data,
                self.INPUT_TABLE_NAME: inputs,
                self.RESULT_TABLE_NAME: results,
            }
        else:
            return {}

    def get_inputs(self, run_id: int, names: list[str] | None = None) -> dict:
        """
        Retrieve input data for a specific run ID from the input table.

        Args:
            run_id (int): The unique identifier for the run whose inputs are to be retrieved.
            names (list[str] | None): A list of input names to filter the results. If None, all inputs are retrieved.

        Returns:
            dict: A dictionary containing the input data, where keys are input names and values are the corresponding data.
        """
        sql_mapping = lambda value: (self._from_json(value))
        return self._get(self.input_table, run_id, names, sql_mapping)

    def get_results(self, run_id: int, names: list[str] | None = None) -> dict:
        """
        Retrieve results from the database for a specific run ID and optional list of names.

        Args:
            run_id (int): The unique identifier for the run whose results are to be retrieved.
            names (list[str] | None): A list of column names to filter the results. If None, all columns are retrieved.

        Returns:
            dict: A dictionary containing the retrieved results, where the keys are column names and the values
                  are the corresponding data, processed using the `_from_bytes` method.
        """
        sql_mapping = lambda value: (self._from_bytes(value))
        return self._get(self.result_table, run_id, names, sql_mapping)

    def delete_run(self, run_id: int) -> None:
        """
        Deletes a run entry from the database based on the provided run ID.

        Args:
            run_id (int): The unique identifier of the run to be deleted.

        Returns:
            None
        """
        self._delete_blob_data_on_disk(run_id)
        condition = self.run_table.c[self.TableKeys.PRIMARY_KEY.value] == run_id
        statement = sql.delete(self.run_table).where(condition)
        self._execute(statement)

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
            self.RUN_TABLE_NAME,
            self.table_metadata,
            sql.Column(
                self.TableKeys.PRIMARY_KEY.value,
                sql.Integer,
                primary_key=True,
                autoincrement=True,
            ),
            sql.Column(
                self.TableKeys.THEORY.value,
                sql.String,
                nullable=False,
            ),
            sql.Column(
                self.TableKeys.COUPLING.value,
                sql.Float,
                nullable=False,
            ),
            sql.Column(
                self.TableKeys.DEGENERACY.value,
                sql.Float,
                nullable=False,
            ),
            sql.Column(
                self.TableKeys.DATE.value,
                sql.String,
                nullable=False,
            ),
            sql.Column(
                self.TableKeys.TIME.value,
                sql.String,
                nullable=False,
            ),
            sql.Column(self.TableKeys.STATUS.value, sql.String, nullable=False),
        )
        self._create_table(table)
        return table

    def _build_inputs_table(self) -> sql.Table:
        """
        Builds and returns the SQLAlchemy Table object for the inputs table.

        This method constructs a table definition for storing input data, using
        the predefined table name and a JSON column type.

        Returns:
            sql.Table: The SQLAlchemy Table object representing the inputs table.
        """
        return self._build_data_table(self.INPUT_TABLE_NAME, sql.JSON)

    def _build_results_table(self) -> sql.Table:
        """
        Constructs and returns the results table for the database.

        This method creates a SQL table with the name specified by
        `RESULTS_TABLE_NAME` and a column of type `LargeBinary` to store
        binary data.

        Returns:
            sql.Table: The constructed results table.
        """
        return self._build_data_table(self.RESULT_TABLE_NAME, sql.LargeBinary)

    def _build_data_table(self, table_name, sql_data_type) -> sql.Table:
        """
        Builds and creates a SQLAlchemy table with the specified name and data type.

        This method defines a table schema with the following columns:
        - `RUN_ID`: An integer column that acts as a foreign key referencing the primary key
          of the runs table. It is non-nullable and enforces cascading deletes.
        - `NAME`: A string column that is non-nullable.
        - `VALUE`: A column with a data type specified by the `sql_data_type` parameter,
          which can be nullable.

        The table also includes a composite primary key constraint on the `RUN_ID` and `NAME` columns.

        After defining the table schema, the method creates the table in the database
        if it does not already exist.

        Args:
            table_name (str): The name of the table to be created.
            sql_data_type (sqlalchemy.types.TypeEngine): The SQLAlchemy data type for the `VALUE` column.

        Returns:
            sqlalchemy.Table: The created SQLAlchemy table object.
        """
        table = sql.Table(
            table_name,
            self.table_metadata,
            sql.Column(
                self.TableKeys.RUN_ID.value,
                sql.Integer,
                sql.ForeignKey(
                    f"{self.RUN_TABLE_NAME}.{self.TableKeys.PRIMARY_KEY.value}",
                    ondelete="CASCADE",
                ),
                nullable=False,
            ),
            sql.Column(
                self.TableKeys.NAME.value,
                sql.String,
                nullable=False,
            ),
            sql.Column(
                self.TableKeys.VALUE.value,
                sql_data_type,
                nullable=True,
            ),
            sql.PrimaryKeyConstraint(
                self.TableKeys.RUN_ID.value, self.TableKeys.NAME.value
            ),
            sql.Index(f"idx_{table_name}_run_id", self.TableKeys.RUN_ID.value),
            sql.Index(f"idx_{table_name}_name", self.TableKeys.NAME.value),
        )
        self._create_table(table)
        return table

    def _create_table(self, table):
        table.create(self.engine, checkfirst=True)

    def _insert_run(self, inputs: any, status: RunStatus):
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
            self.TableKeys.THEORY.value: inputs.theory,
            self.TableKeys.COUPLING.value: inputs.coupling,
            self.TableKeys.DEGENERACY.value: inputs.degeneracy,
            self.TableKeys.DATE.value: now.date().isoformat(),
            self.TableKeys.TIME.value: now.time().isoformat(),
            self.TableKeys.STATUS.value: status.value,
        }
        statement = sql.insert(self.run_table).values(data)
        result = self._execute(statement)
        if run_id := result.inserted_primary_key:
            self.run_id = run_id[0]

    def _delete_blob_data_on_disk(self, run_id: int):
        try:
            fixed_table = sql.Table(
                self.FIXED_TABLE_NAME, self.table_metadata, autoload_with=self.engine
            )
        except NoSuchTableError:
            return
        stmt_select = sql.select(fixed_table.c.adr).where(
            fixed_table.c.run_id == run_id
        )
        result = self._execute(stmt_select).fetchall()
        file_paths = [row[0] for row in result]
        for path in file_paths:
            try:
                os.remove(path)
            except FileNotFoundError:
                print(f"Warning: file not found: {path}")
            except Exception as e:
                print(f"Warning: failed to delete {path}: {e}")

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

    def _insert_from_dict(
        self,
        table,
        data: dict[str, any],
        sql_mapping: Callable[[any], any],
        conflict_mode: ConflictMode = ConflictMode.FAIL,
    ) -> None:
        """
        Inserts data into a specified table by mapping values through a provided SQL mapping function.

        Args:
            table (str): The name of the table where the data will be inserted.
            data (dict[str, any]): A dictionary containing column names as keys and their corresponding values.
            sql_mapping (Callable[[any], any]): A function that maps input values to their SQL-compatible representations.

        Returns:
            None
        """
        for name, value in data.items():
            if mapped_value := sql_mapping(value):
                self._insert(table, name, mapped_value, conflict_mode)

    def _insert(
        self,
        table: sql.Table,
        name: str,
        value: any,
        conflict_mode: ConflictMode = ConflictMode.FAIL,
    ):
        """
        Inserts a record into the specified SQL table with the given name and value, handling conflicts according to the specified mode.
        Args:
            table (sql.Table): The SQLAlchemy table object where the record will be inserted.
            name (str): The name/key associated with the value to insert.
            value (any): The value to be inserted into the table.
            conflict_mode (ConflictMode, optional): Specifies how to handle conflicts on unique constraints.
                Defaults to ConflictMode.FAIL. If set to ConflictMode.UPDATE, existing records with the same
                run_id and name will be updated with the new value.
        Returns:
            None
        Raises:
            Any exceptions raised by the underlying database execution.
        """
        data = {
            self.TableKeys.RUN_ID.value: self.run_id,
            self.TableKeys.NAME.value: name,
            self.TableKeys.VALUE.value: value,
        }
        statement = sqlite_insert(table).values(data)
        if conflict_mode == self.ConflictMode.UPDATE:
            statement = statement.on_conflict_do_update(
                index_elements=[
                    self.TableKeys.RUN_ID.value,
                    self.TableKeys.NAME.value,
                ],
                set_={self.TableKeys.VALUE.value: value},
            )
        self._execute(statement)

    def _insert_with_update(self, table: sql.Table, name: str, value: any):
        """
        Inserts a record into the specified SQL table or updates it if a conflict occurs.

        Args:
            table (sql.Table): The SQLAlchemy Table object representing the target table.
            name (str): The name of the record to insert or update.
            value (any): The value associated with the record.

        Behavior:
            - If a record with the same `RUN_ID` and `NAME` already exists in the table,
              the `VALUE` field of the existing record will be updated.
            - If no such record exists, a new record will be inserted.

        Raises:
            Any exceptions raised by the `_execute` method during the execution of the SQL statement.
        """
        data = {
            self.TableKeys.RUN_ID.value: self.run_id,
            self.TableKeys.NAME.value: name,
            self.TableKeys.VALUE.value: value,
        }
        statement = (
            sqlite_insert(table)
            .values(data)
            .on_conflict_do_update(
                index_elements=[
                    self.TableKeys.RUN_ID.value,
                    self.TableKeys.NAME.value,
                ],
                set_={self.TableKeys.VALUE.value: value},
            )
        )
        self._execute(statement)

    def _get(
        self,
        table: sql.Table,
        run_id: int,
        names: list[str] | None,
        sql_mapping: Callable[[any], any],
    ) -> dict:
        """
        Retrieve data from a specified SQL table based on a run ID and optional list of names.

        Args:
            table (sql.Table): The SQLAlchemy Table object to query.
            run_id (int): The run ID to filter the data.
            names (list[str] | None): An optional list of names to filter the data. If None, no name filtering is applied.
            sql_mapping (Callable[[any], any]): A callable to transform the SQL value into the desired format.

        Returns:
            dict: A dictionary where the keys are the names from the table and the values are the transformed values
                  obtained by applying `sql_mapping` to the corresponding SQL values.

        """
        conditions = [table.c[self.TableKeys.RUN_ID.value] == run_id]
        if names is not None:
            conditions.append(table.c[self.TableKeys.NAME.value].in_(names))
        statement = sql.select(table).where(*conditions)
        rows = self._execute(statement).mappings().all()
        return {
            row[self.TableKeys.NAME.value]: sql_mapping(row[self.TableKeys.VALUE.value])
            for row in rows
        }

    def _execute(self, statement) -> sql.CursorResult[any]:
        """
        Executes a given SQL statement using the database engine.

        This method establishes a connection to the database, executes the provided
        SQL statement, and returns the result.

        Args:
            statement: The SQL statement to be executed.

        Returns:
            sql.CursorResult[any]: The result of the executed SQL statement.
        """
        with self.engine.begin() as connection:
            result = connection.execute(statement)
            return result

    def _to_bytes(self, data: float | np.ndarray) -> bytes | None:
        """
        Converts a float or a NumPy array into a bytes representation.

        Parameters:
            data (float | np.ndarray): The input data to be converted. It can be either
                a float or a NumPy array.

        Returns:
            bytes | None: The bytes representation of the input data if it is a float
                or a NumPy array. Returns None if the input data type is unsupported.
        """
        if isinstance(data, float):
            return struct.pack("d", data)
        elif isinstance(data, np.ndarray):
            arr_bytes = io.BytesIO()
            np.save(arr_bytes, data)
            compressed_arr_bytes = blosc2.compress(arr_bytes.getvalue())
            return compressed_arr_bytes
        else:
            return None

    def _from_bytes(self, data: bytes) -> float | np.ndarray | None:
        """
        Converts a byte sequence into a float, a NumPy array, or None.

        This method attempts to interpret the input byte sequence as either:
        - A double-precision floating-point number if the length of the data is 8 bytes.
        - A NumPy array if the data represents a serialized array.
        - Returns None if the conversion fails.

        Args:
            data (bytes): The byte sequence to be converted.

        Returns:
            float | np.ndarray | None: The converted value as a float, a NumPy array,
            or None if the conversion is unsuccessful.
        """
        try:
            if len(data) == 8:
                return struct.unpack("d", data)[0]
            else:
                decompressed_data = blosc2.decompress(data)
                arr_bytes = io.BytesIO(decompressed_data)
                return np.load(arr_bytes, allow_pickle=False)
        except Exception:
            return None

    def _to_json(self, data: any) -> json:
        """
        Converts the given data to a JSON-formatted string.

        Args:
            data (any): The data to be converted to JSON.

        Returns:
            str: A JSON-formatted string representation of the data if conversion is successful.
            None: If an error occurs during the conversion process.
        """
        try:
            return json.dumps(data)
        except:
            return None

    def _from_json(self, data: json) -> any:
        """
        Converts a JSON-formatted string into a Python object.

        Args:
            data (json): A JSON-formatted string to be deserialized.

        Returns:
            any: The deserialized Python object if the input is valid JSON.
                 Returns None if deserialization fails.
        """
        try:
            return json.loads(data)
        except:
            return None
