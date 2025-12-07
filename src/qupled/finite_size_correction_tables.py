from datetime import datetime
import sqlalchemy as sql

from . import native
from . import base_tables as base

INPUT_TABLE_NAME = "fsc_inputs"
RESULT_TABLE_NAME = "fsc_results"
RUN_TABLE_NAME = "fsc_runs"


class FiniteSizeCorrectionTables(base.BaseTables):
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
                base.TableKeys.NUMBER_OF_PARTICLES.value,
                sql.Integer,
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

    def _insert_run(self, inputs: any, status: base.RunStatus):
        """ """
        inputs_scheme = inputs.scheme
        now = datetime.now()
        data = {
            base.TableKeys.THEORY.value: inputs_scheme.theory,
            base.TableKeys.COUPLING.value: inputs_scheme.coupling,
            base.TableKeys.DEGENERACY.value: inputs_scheme.degeneracy,
            base.TableKeys.NUMBER_OF_PARTICLES.value: inputs.number_of_particles,
            base.TableKeys.DATE.value: now.date().isoformat(),
            base.TableKeys.TIME.value: now.time().isoformat(),
            base.TableKeys.STATUS.value: status.value,
        }
        statement = sql.insert(self.run_table).values(data)
        result = self._execute(statement)
        if run_id := result.inserted_primary_key:
            self.run_id = run_id[0]
