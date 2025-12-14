import io
import json
import struct
from enum import Enum
from collections.abc import Callable

import numpy as np
import sqlalchemy as sql
from sqlalchemy.dialects.sqlite import insert as sqlite_insert
import blosc2


class TableKeys(Enum):
    NAME = "name"
    PRIMARY_KEY = "id"
    RUN_ID = "run_id"
    STATUS = "status"
    VALUE = "value"


class RunStatus(Enum):
    RUNNING = "STARTED"
    SUCCESS = "SUCCESS"
    FAILED = "FAILED"


class ConflictMode(Enum):
    FAIL = "FAIL"
    UPDATE = "UPDATE"


class BaseTables:
    """
    Base class for managing database tables, including run, input, and result tables.
    Provides methods for inserting, retrieving, and updating data in the database.
    """

    def __init__(
        self,
        engine: sql.Engine,
        run_table_name: str,
        input_table_name: str,
        result_table_name: str,
    ):
        """
        Initializes the BaseTables instance with the database engine and table names.

        Args:
            engine (sql.Engine): SQLAlchemy engine for database connection.
            run_table_name (str): Name of the run table.
            input_table_name (str): Name of the input table.
            result_table_name (str): Name of the result table.
        """
        # Set sqlite properties
        self.engine = engine
        # Store table names
        self.run_table_name = run_table_name
        self.input_table_name = input_table_name
        self.result_table_name = result_table_name
        # Create tables
        self.table_metadata = sql.MetaData()
        self.input_table: sql.Table | None = None
        self.result_table: sql.Table | None = None
        self.run_table: sql.Table | None = None
        self.run_id: int | None = None

    def insert_run(self, inputs):
        """
        Inserts a new run into the database and adds its associated inputs.

        Args:
            inputs: Input data for the run.
        """
        self._insert_run(inputs, RunStatus.RUNNING)
        self.insert_inputs(inputs.__dict__)

    def insert_inputs(self, inputs: dict[str, any]):
        """
        Inserts input data into the input table for the current run.

        Args:
            inputs (dict[str, any]): Dictionary of input data to insert.
        """
        if self.run_id is not None:
            sql_mapping = lambda value: (self._to_json(value))
            self._insert_data_from_dict(self.input_table, inputs, sql_mapping)

    def insert_results(
        self,
        results: dict[str, any],
        conflict_mode: ConflictMode = ConflictMode.FAIL,
    ):
        """
        Inserts result data into the result table for the current run.

        Args:
            results (dict[str, any]): Dictionary of result data to insert.
            conflict_mode (ConflictMode): Conflict resolution mode (default is FAIL).
        """
        if self.run_id is not None:
            sql_mapping = lambda value: (self._to_bytes(value))
            self._insert_data_from_dict(
                self.result_table, results, sql_mapping, conflict_mode
            )

    def inspect_runs(self) -> list[dict[str, any]]:
        """
        Retrieves all runs from the run table.

        Returns:
            list[dict[str, any]]: List of dictionaries containing run data.
        """
        statement = sql.select(self.run_table)
        rows = self._execute(statement).mappings().all()
        return [{key: row[key] for key in row.keys()} for row in rows]

    def update_run_status(self, status: RunStatus) -> None:
        """
        Updates the status of the current run in the run table.

        Args:
            status (RunStatus): New status to set for the run.
        """
        if self.run_id is not None:
            statement = (
                sql.update(self.run_table)
                .where(self.run_table.c[TableKeys.PRIMARY_KEY.value] == self.run_id)
                .values({TableKeys.STATUS.value: status.value})
            )
            self._execute(statement)

    def get_run(
        self,
        run_id: int,
        input_names: list[str] | None = None,
        result_names: list[str] | None = None,
    ) -> dict:
        """
        Retrieves data for a specific run, including its inputs and results.

        Args:
            run_id (int): ID of the run to retrieve.
            input_names (list[str] | None): Optional list of input names to filter.
            result_names (list[str] | None): Optional list of result names to filter.

        Returns:
            dict: Dictionary containing run, input, and result data.
        """
        statement = sql.select(self.run_table).where(
            self.run_table.c[TableKeys.PRIMARY_KEY.value] == run_id
        )
        result = self._execute(statement).mappings().first()
        if result is not None:
            run_data = {key: result[key] for key in result.keys()}
            inputs = self.get_inputs(run_id, names=input_names)
            results = self.get_results(run_id, names=result_names)
            return {
                self.run_table_name: run_data,
                self.input_table_name: inputs,
                self.result_table_name: results,
            }
        return {}

    def get_inputs(self, run_id: int, names: list[str] | None = None) -> dict:
        """
        Retrieves input data for a specific run.

        Args:
            run_id (int): ID of the run to retrieve inputs for.
            names (list[str] | None): Optional list of input names to filter.

        Returns:
            dict: Dictionary of input data.
        """
        sql_mapping = lambda value: (self._from_json(value))
        return self._get_data(self.input_table, run_id, names, sql_mapping)

    def get_results(self, run_id: int, names: list[str] | None = None) -> dict:
        """
        Retrieves result data for a specific run.

        Args:
            run_id (int): ID of the run to retrieve results for.
            names (list[str] | None): Optional list of result names to filter.

        Returns:
            dict: Dictionary of result data.
        """
        sql_mapping = lambda value: (self._from_bytes(value))
        return self._get_data(self.result_table, run_id, names, sql_mapping)

    def delete_run(self, run_id: int) -> None:
        """
        Deletes a specific run from the database.

        Args:
            run_id (int): ID of the run to delete.
        """
        condition = self.run_table.c[TableKeys.PRIMARY_KEY.value] == run_id
        statement = sql.delete(self.run_table).where(condition)
        self._execute(statement)

    def _build_tables(self):
        """
        Builds the run, input, and result tables in the database.
        """
        self.run_table = self._build_run_table()
        self.input_table = self._build_inputs_table()
        self.result_table = self._build_results_table()

    def _build_run_table(self) -> sql.Table:
        """
        Builds the run table.

        Returns:
            sql.Table: SQLAlchemy Table object for the run table.

        Raises:
            NotImplementedError: If not implemented in a subclass.
        """
        raise NotImplementedError("This method should be implemented in a subclass.")

    def _build_inputs_table(self) -> sql.Table:
        """
        Builds the input table.

        Returns:
            sql.Table: SQLAlchemy Table object for the input table.
        """
        return self._build_data_table(self.input_table_name, sql.JSON)

    def _build_results_table(self) -> sql.Table:
        """
        Builds the result table.

        Returns:
            sql.Table: SQLAlchemy Table object for the result table.
        """
        return self._build_data_table(self.result_table_name, sql.LargeBinary)

    def _build_data_table(self, table_name, sql_data_type) -> sql.Table:
        """
        Builds a generic data table with the specified name and data type.

        Args:
            table_name (str): Name of the table to create.
            sql_data_type: SQLAlchemy data type for the table's value column.

        Returns:
            sql.Table: SQLAlchemy Table object for the data table.
        """
        table = sql.Table(
            table_name,
            self.table_metadata,
            sql.Column(
                TableKeys.RUN_ID.value,
                sql.Integer,
                sql.ForeignKey(
                    f"{self.run_table_name}.{TableKeys.PRIMARY_KEY.value}",
                    ondelete="CASCADE",
                ),
                nullable=False,
            ),
            sql.Column(
                TableKeys.NAME.value,
                sql.String,
                nullable=False,
            ),
            sql.Column(
                TableKeys.VALUE.value,
                sql_data_type,
                nullable=True,
            ),
            sql.PrimaryKeyConstraint(TableKeys.RUN_ID.value, TableKeys.NAME.value),
            sql.Index(f"idx_{table_name}_run_id", TableKeys.RUN_ID.value),
            sql.Index(f"idx_{table_name}_name", TableKeys.NAME.value),
        )
        self._create_table(table)
        return table

    def _create_table(self, table: sql.Table):
        """
        Creates a table in the database if it does not already exist.

        Args:
            table (sql.Table): SQLAlchemy Table object to create.
        """
        table.create(self.engine, checkfirst=True)

    def _insert_run(self, inputs: any, status: RunStatus):
        """
        Inserts a new run into the run table.

        Args:
            inputs (any): Input data for the run.
            status (RunStatus): Status of the run.

        Raises:
            NotImplementedError: If not implemented in a subclass.
        """
        raise NotImplementedError("This method should be implemented in a subclass.")

    def _insert_data_from_dict(
        self,
        table,
        data: dict[str, any],
        sql_mapping: Callable[[any], any],
        conflict_mode: ConflictMode = ConflictMode.FAIL,
    ) -> None:
        """
        Inserts data into a table from a dictionary.

        Args:
            table: SQLAlchemy Table object to insert data into.
            data (dict[str, any]): Dictionary of data to insert.
            sql_mapping (Callable[[any], any]): Mapping function to transform values.
            conflict_mode (ConflictMode): Conflict resolution mode (default is FAIL).
        """
        for name, value in data.items():
            if mapped_value := sql_mapping(value):
                self._insert_data(table, name, mapped_value, conflict_mode)

    def _insert_data(
        self,
        table: sql.Table,
        name: str,
        value: any,
        conflict_mode: ConflictMode = ConflictMode.FAIL,
    ):
        """
        Inserts a single data entry into a table.

        Args:
            table (sql.Table): SQLAlchemy Table object to insert data into.
            name (str): Name of the data entry.
            value (any): Value of the data entry.
            conflict_mode (ConflictMode): Conflict resolution mode (default is FAIL).
        """
        data = {
            TableKeys.RUN_ID.value: self.run_id,
            TableKeys.NAME.value: name,
            TableKeys.VALUE.value: value,
        }
        statement = sqlite_insert(table).values(data)
        if conflict_mode == ConflictMode.UPDATE:
            statement = statement.on_conflict_do_update(
                index_elements=[
                    TableKeys.RUN_ID.value,
                    TableKeys.NAME.value,
                ],
                set_={TableKeys.VALUE.value: value},
            )
        self._execute(statement)

    def _get_data(
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
        conditions = [table.c[TableKeys.RUN_ID.value] == run_id]
        if names is not None:
            conditions.append(table.c[TableKeys.NAME.value].in_(names))
        statement = sql.select(table).where(*conditions)
        rows = self._execute(statement).mappings().all()
        return {
            row[TableKeys.NAME.value]: sql_mapping(row[TableKeys.VALUE.value])
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
            if hasattr(data, "to_dict") and callable(data.to_dict):
                return json.dumps(data.to_dict())
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
