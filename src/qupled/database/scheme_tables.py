from datetime import datetime
from enum import Enum

import sqlalchemy as sql

from qupled.database.base_tables import (
    BaseTables,
    RunStatus,
    TableKeys as BaseTableKeys,
)
from qupled.native import delete_blob_data_on_disk

INPUT_TABLE_NAME = "inputs"
RESULT_TABLE_NAME = "results"
RUN_TABLE_NAME = "runs"


class TableKeys(Enum):
    COUPLING = "coupling"
    DATE = "date"
    DEGENERACY = "degeneracy"
    THEORY = "theory"
    TIME = "time"


class SchemeTables(BaseTables):
    """ """

    def __init__(self, engine: sql.Engine):
        """ """
        super().__init__(engine, RUN_TABLE_NAME, INPUT_TABLE_NAME, RESULT_TABLE_NAME)
        self._build_tables()

    def delete_run(self, run_id: int) -> None:
        """ """
        delete_blob_data_on_disk(self.engine.url.database, run_id)
        super().delete_run(run_id)

    def _build_run_table(self):
        """ """
        table = sql.Table(
            self.run_table_name,
            self.table_metadata,
            sql.Column(
                BaseTableKeys.PRIMARY_KEY.value,
                sql.Integer,
                primary_key=True,
                autoincrement=True,
            ),
            sql.Column(
                TableKeys.THEORY.value,
                sql.String,
                nullable=False,
            ),
            sql.Column(
                TableKeys.COUPLING.value,
                sql.Float,
                nullable=False,
            ),
            sql.Column(
                TableKeys.DEGENERACY.value,
                sql.Float,
                nullable=False,
            ),
            sql.Column(
                TableKeys.DATE.value,
                sql.String,
                nullable=False,
            ),
            sql.Column(
                TableKeys.TIME.value,
                sql.String,
                nullable=False,
            ),
            sql.Column(
                BaseTableKeys.STATUS.value,
                sql.String,
                nullable=False,
            ),
        )
        self._create_table(table)
        return table

    def _insert_run(self, inputs: any, status: RunStatus):
        """ """
        now = datetime.now()
        data = {
            TableKeys.THEORY.value: inputs.theory,
            TableKeys.COUPLING.value: inputs.coupling,
            TableKeys.DEGENERACY.value: inputs.degeneracy,
            TableKeys.DATE.value: now.date().isoformat(),
            TableKeys.TIME.value: now.time().isoformat(),
            BaseTableKeys.STATUS.value: status.value,
        }
        statement = sql.insert(self.run_table).values(data)
        result = self._execute(statement)
        if run_id := result.inserted_primary_key:
            self.run_id = run_id[0]
