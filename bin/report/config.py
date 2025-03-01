"""Configuration of the report module."""

import os
from datetime import datetime
from pathlib import Path

GITHUB_URL = 'https://github.com/maelyg/ont_amplicon'


class Config:

    METADATA_FILE = 'index.csv'
    TIMESTAMP_FILE = 'timestamp.txt'
    REPORT_FILE = 'report.html'

    class REPORT:
        TITLE = "Amplicon sequencing assembly report"
        SUBTITLE = (
            "Results generated through the"
            f' <a href="{GITHUB_URL}" target="_blank">'
            'ONT-amplicon NextFlow workflow</a>.')

    class CRITERIA:
        MIN_RAW_READS = 4000
        MIN_FILTERED_READS = 1000

    def load_results(self, result_dir: Path):
        """Read the results from the result directory."""
        os.environ['RESULT_DIR'] = str(result_dir)
        self.sample_id = self._get_sample_id(result_dir)
        self.start_time = self._read_timestamp(result_dir)

    @property
    def result_dir(self) -> Path:
        """Ensure that result dir is propagated between Config instances."""
        return Path(os.environ['RESULT_DIR'])

    @property
    def report_path(self) -> Path:
        return self.result_dir / self.REPORT_FILE

    @property
    def metadata_path(self) -> Path:
        return self.result_dir / self.METADATA_FILE

    @property
    def parameters_path(self) -> Path:
        return self.result_dir / self.PARAMS_FILE  # TODO: not available yet

    @property
    def run_stats_path(self) -> Path:
        return list(self.result_dir.glob('run_qc_report_*.txt'))[0]

    @property
    def run_qc_html_file(self) -> Path:
        return list(self.result_dir.glob('run_qc_report_*.html'))[0].name

    def _get_sample_id(self, result_dir: Path):
        """Return the sample ID from the result directory."""
        bam_path = list(result_dir.glob('*.bam'))[0]
        return bam_path.name.split('_aln.')[0]

    def _read_timestamp(self, result_dir: Path) -> str:
        """Return the timestamp of the start of the workflow."""
        timestamp_path = result_dir / self.TIMESTAMP_FILE
        if timestamp_path.exists():
            return datetime.strptime(
                timestamp_path.read_text().strip(),
                '%Y-%m-%d %H:%M:%S',
            )
        return None
