# -----------------------------------------------------------------------
# DataBaseHandler class
# -----------------------------------------------------------------------

import io
import json
import struct
from datetime import datetime
from enum import Enum
from typing import Any, Callable

import numpy as np
import sqlalchemy as sql


class DataBaseHandler:

    DATABASE_NAME = "scheme_results.db"
    RUNS_TABLE_NAME = "Runs"
    INPUTS_TABLE_NAME = "Inputs"
    RESULTS_TABLE_NAME = "Results"

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

    def __init__(self):
        self.engine = sql.create_engine(f"sqlite:///{self.DATABASE_NAME}")
        self.table_metadata = sql.MetaData()
        self.run_table = self._build_run_table()
        self.inputs_table = self._build_inputs_table()
        self.results_table = self._build_results_table()
        self.run_id: int | None = None

    def insert_run(self, inputs, results):
        self.insert_run_data(inputs)
        self.insert_inputs_data(inputs)
        self.insert_results_data(results)
        return self.run_id

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
                self.TableKeys.PRIMARY_KEY.value,
                sql.Integer,
                primary_key=True,
                autoincrement=True,
            ),
            sql.Column(
                self.TableKeys.RUN_ID.value,
                sql.Integer,
                sql.ForeignKey(
                    f"{self.RUNS_TABLE_NAME}.{self.TableKeys.PRIMARY_KEY.value}"
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
        )
        table.create(self.engine, checkfirst=True)
        return table

    def insert_run_data(self, inputs: Any):
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

    def insert_inputs_data(self, inputs: Any):
        sql_mapping = lambda value: (self._to_json(value))
        self._insert_data(self.inputs_table, inputs, sql_mapping)

    def insert_results_data(self, results: Any):
        sql_mapping = lambda value: (self._to_bytes(value))
        self._insert_data(self.results_table, results, sql_mapping)

    def _insert_data(
        self, table: sql.Table, data_src: Any, sql_mapping: Callable[[Any], Any]
    ):
        for attr, value in data_src.__dict__.items():
            if mapped_value := sql_mapping(value):
                data = {
                    self.TableKeys.RUN_ID.value: self.run_id,
                    self.TableKeys.NAME.value: attr,
                    self.TableKeys.VALUE.value: mapped_value,
                }
                statement = sql.insert(table).values(data)
                self._execute(statement)

    def _execute(self, statement) -> sql.CursorResult[Any]:
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

    def _to_json(self, data: any) -> json:
        try:
            return json.dumps(data)
        except:
            return None
