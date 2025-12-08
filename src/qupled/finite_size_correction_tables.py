from datetime import datetime
from enum import Enum

import sqlalchemy as sql

from .base_tables import BaseTables, RunStatus, TableKeys as BaseTableKeys

INPUT_TABLE_NAME = "fsc_inputs"
RESULT_TABLE_NAME = "fsc_results"
RUN_TABLE_NAME = "fsc_runs"


class TableKeys(Enum):
    COUPLING = "coupling"
    DATE = "date"
    DEGENERACY = "degeneracy"
    NUMBER_OF_PARTICLES = "number_of_particles"
    SCHEME_RUN_ID = "scheme_run_id"
    THEORY = "theory"
    TIME = "time"


class FiniteSizeCorrectionTables(BaseTables):
    """ """

    def __init__(self, engine: sql.Engine):
        """ """
        super().__init__(engine, RUN_TABLE_NAME, INPUT_TABLE_NAME, RESULT_TABLE_NAME)
        self._build_tables()

    def insert_inputs(self, inputs: dict[str, any]):
        """ """
        inputs_local = dict(inputs)
        inputs_local.pop("scheme", None)
        super().insert_inputs(inputs_local)

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
                TableKeys.NUMBER_OF_PARTICLES.value,
                sql.Integer,
                nullable=False,
            ),
            sql.Column(
                TableKeys.SCHEME_RUN_ID.value,
                sql.String,
                nullable=True,
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

    def update_scheme_run_id(self, scheme_run_id: int) -> None:
        """ """
        if self.run_id is not None:
            statement = (
                sql.update(self.run_table)
                .where(self.run_table.c[BaseTableKeys.PRIMARY_KEY.value] == self.run_id)
                .values({TableKeys.SCHEME_RUN_ID.value: scheme_run_id})
            )
            self._execute(statement)

    def _insert_run(self, inputs: any, status: RunStatus):
        """ """
        inputs_scheme = inputs.scheme
        now = datetime.now()
        data = {
            TableKeys.THEORY.value: inputs_scheme.theory,
            TableKeys.COUPLING.value: inputs_scheme.coupling,
            TableKeys.DEGENERACY.value: inputs_scheme.degeneracy,
            TableKeys.NUMBER_OF_PARTICLES.value: inputs.number_of_particles,
            TableKeys.DATE.value: now.date().isoformat(),
            TableKeys.TIME.value: now.time().isoformat(),
            BaseTableKeys.STATUS.value: status.value,
        }
        statement = sql.insert(self.run_table).values(data)
        result = self._execute(statement)
        if run_id := result.inserted_primary_key:
            self.run_id = run_id[0]
