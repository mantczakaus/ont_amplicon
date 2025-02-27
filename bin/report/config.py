"""Configuration of the report module."""

from pathlib import Path

GITHUB_URL = 'https://github.com/maelyg/ont_amplicon'


class Config:

    class REPORT:
        TITLE = "Amplicon sequencing assembly report"
        SUBTITLE = (
            "Results generated through the"
            f' <a href="{GITHUB_URL}" target="_blank">'
            'ONT-amplicon NextFlow workflow</a>.')

    @property
    def metadata_path(self) -> Path:
        return self.result_dir / 'index.csv'

    def read_results(self, result_dir: Path):
        """Read the results from the result directory."""
        self.start_time = self._read_timestamp(result_dir)
        self.sample_id = self._get_sample_id(result_dir)

    def _get_sample_id(self, result_dir: Path):
        """Return the sample ID from the result directory."""
        bam_path = list(result_dir.glob('*.bam'))[0]
        return bam_path.split('_aln.')[0].name

    def _read_timestamp(result_dir: Path) -> str:
        """Return the timestamp of the start of the workflow."""
        timestamp_path = result_dir / 'timestamp.txt'
        if timestamp_path.exists():
            return timestamp_path.read_text().strip()
        return None
