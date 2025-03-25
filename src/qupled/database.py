# -----------------------------------------------------------------------
# DataBaseHandler class
# -----------------------------------------------------------------------

import io
import json
import struct
from datetime import datetime
from enum import Enum
from collections.abc import Callable

import numpy as np
import sqlalchemy as sql
from sqlalchemy.dialects.sqlite import insert as sqlite_insert


class DataBaseHandler:

    DEFAULT_DATABASE_NAME = "scheme_results.db"
    RUNS_TABLE_NAME = "runs"
    INPUTS_TABLE_NAME = "inputs"
    RESULTS_TABLE_NAME = "results"

    class TableKeys(Enum):
        COUPLING = "coupling"
        DATE = "date"
        DEGENERACY = "degeneracy"
        NAME = "name"
        PRIMARY_KEY = "id"
        RUN_ID = "run_id"
        THEORY = "theory"
        TIME = "time"
        VALUE = "value"

    def __init__(self, database_name: str | None = None):
        self.database_name = (
            database_name if database_name is not None else self.DEFAULT_DATABASE_NAME
        )
        self.engine = sql.create_engine(f"sqlite:///{self.database_name}")
        DataBaseHandler._set_sqlite_pragma(self.engine) # Enforce foreign keys in sqlite
        self.table_metadata = sql.MetaData()
        self.run_table = self._build_run_table()
        self.input_table = self._build_inputs_table()
        self.result_table = self._build_results_table()
        self.run_id: int | None = None

    def insert_run(self, inputs, results):
        self._insert_run(inputs)
        self.insert_inputs(inputs.__dict__)
        self.insert_results(results.__dict__)

    def insert_inputs(self, inputs: dict[str, any]):
        sql_mapping = lambda value: (self._to_json(value))
        self._insert_from_dict(self.input_table, inputs, sql_mapping)

    def insert_results(self, results: dict[str, any]):
        sql_mapping = lambda value: (self._to_bytes(value))
        self._insert_from_dict(self.result_table, results, sql_mapping)

    def inspect_runs(self) -> list[dict[str, any]]:
        statement = sql.select(self.run_table)
        rows = self._execute(statement).mappings().all()
        return [{key: row[key] for key in row.keys()} for row in rows]

    def get_run(self, run_id: int) -> dict:
        statement = sql.select(self.run_table).where(
            self.run_table.c[self.TableKeys.PRIMARY_KEY.value] == run_id
        )
        result = self._execute(statement).mappings().first()
        if result is not None: 
            run_data = {key: result[key] for key in result.keys()}
            inputs = self.get_inputs(run_id, names=None)
            results = self.get_results(run_id, names=None)
            return {
                self.RUNS_TABLE_NAME: run_data,
                self.INPUTS_TABLE_NAME: inputs,
                self.RESULTS_TABLE_NAME: results
            }
        else: 
            return {}

    def get_inputs(self, run_id: int, names: list[str] | None) -> dict:
        sql_mapping = lambda value: (self._from_json(value))
        return self._get(self.input_table, run_id, names, sql_mapping)

    def get_results(self, run_id: int, names: list[str] | None) -> dict:
        sql_mapping = lambda value: (self._from_bytes(value))
        return self._get(self.result_table, run_id, names, sql_mapping)

    def delete_run(self, run_id: int) -> None:
        condition = self.run_table.c[self.TableKeys.PRIMARY_KEY.value] == run_id
        statement = sql.delete(self.run_table).where(condition)
        self._execute(statement)

    def _build_run_table(self):
        table = sql.Table(
            self.RUNS_TABLE_NAME,
            self.table_metadata,
            sql.Column(
                self.TableKeys.PRIMARY_KEY.value,
                sql.Integer,
                primary_key=True,
                autoincrement=True,
            ),
            sql.Column(
                self.TableKeys.THEORY.value,
                sql.JSON,
                nullable=False,
            ),
            sql.Column(
                self.TableKeys.COUPLING.value,
                sql.JSON,
                nullable=False,
            ),
            sql.Column(
                self.TableKeys.DEGENERACY.value,
                sql.JSON,
                nullable=False,
            ),
            sql.Column(
                self.TableKeys.DATE.value,
                sql.JSON,
                nullable=False,
            ),
            sql.Column(
                self.TableKeys.TIME.value,
                sql.JSON,
                nullable=False,
            ),
        )
        table.create(self.engine, checkfirst=True)
        return table

    def _build_inputs_table(self) -> sql.Table:
        return self._build_data_table(self.INPUTS_TABLE_NAME, sql.JSON)

    def _build_results_table(self) -> sql.Table:
        return self._build_data_table(self.RESULTS_TABLE_NAME, sql.LargeBinary)

    def _build_data_table(self, table_name, sql_data_type) -> sql.Table:
        table = sql.Table(
            table_name,
            self.table_metadata,
            sql.Column(
                self.TableKeys.RUN_ID.value,
                sql.Integer,
                sql.ForeignKey(
                    f"{self.RUNS_TABLE_NAME}.{self.TableKeys.PRIMARY_KEY.value}",
                    ondelete="CASCADE"
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
        )
        table.create(self.engine, checkfirst=True)
        return table

    def _insert_run(self, inputs: any):
        now = datetime.now()
        data = {
            self.TableKeys.THEORY.value: self._to_json(inputs.theory),
            self.TableKeys.COUPLING.value: self._to_json(inputs.coupling),
            self.TableKeys.DEGENERACY.value: self._to_json(inputs.degeneracy),
            self.TableKeys.DATE.value: now.date().isoformat(),
            self.TableKeys.TIME.value: now.time().isoformat(),
        }
        statement = sql.insert(self.run_table).values(data)
        result = self._execute(statement)
        if run_id := result.inserted_primary_key:
            self.run_id = run_id[0]

    @staticmethod
    def _set_sqlite_pragma(engine):
        @sql.event.listens_for(engine, "connect")
        def _set_pragma(dbapi_connection, connection_record):
            cursor = dbapi_connection.cursor()
            cursor.execute("PRAGMA foreign_keys=ON")
            cursor.close()

    def _insert_from_dict(
        self, table, data: dict[str, any], sql_mapping: Callable[[any], any]
    ) -> None:
        for name, value in data.items():
            if mapped_value := sql_mapping(value):
                self._insert(table, name, mapped_value)

    def _insert(self, table: sql.Table, name: str, value: any):
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
        with self.engine.begin() as connection:
            result = connection.execute(statement)
            return result

    def _to_bytes(self, data: float | np.ndarray) -> bytes | None:
        if isinstance(data, float):
            return struct.pack("d", data)
        elif isinstance(data, np.ndarray):
            arr_bytes = io.BytesIO()
            np.save(arr_bytes, data)
            return arr_bytes.getvalue()
        else:
            return None

    def _from_bytes(self, data: bytes) -> float | np.ndarray | None:
        try:
            if len(data) == 8:
                return struct.unpack("d", data)[0]
            else:
                arr_bytes = io.BytesIO(data)
                return np.load(arr_bytes, allow_pickle=False)
        except Exception:
            return None

    def _to_json(self, data: any) -> json:
        try:
            return json.dumps(data)
        except:
            return None

    def _from_json(self, data: json) -> any:
        try:
            return json.loads(data)
        except:
            return None
