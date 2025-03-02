"""Define specific results from the analysis."""

import base64

from .config import Config

config = Config()


class FLAGS:
    SUCCESS = 'success'
    WARNING = 'warning'
    DANGER = 'danger'
    NONE = 'secondary'


class AbstractDataRow:

    COLUMNS = []

    def __init__(self, row):
        for colname, _type in self.COLUMNS:
            try:
                value = _type(row[colname].strip())
            except KeyError:
                raise KeyError(f"Column '{colname}' not found in"
                               f" {self.__class__.__name__} row - got:"
                               f" {row.keys()}")
            setattr(self, colname, value)

    def to_json(self):
        return {
            colname: getattr(self, colname)
            for colname, _ in self.COLUMNS
        }


class Metadata(AbstractDataRow):
    """Report the sample metadata."""

    COLUMNS = [
        ('sample_files', str),
        ('target_size', int),
        ('spp_targets', str),
        ('gene_targets', str),
    ]


class RunQC(AbstractDataRow):
    """Report the quality control outcomes."""

    COLUMNS = [
        ('raw_reads', int),
        ('quality_filtered_reads', int),
        ('percent_quality_filtered', float),
        ('raw_reads_flag', str),
        ('qfiltered_flag', str),
    ]

    @property
    def flag(self):
        raw_threshold = self.raw_reads > config.CRITERIA.MIN_RAW_READS
        qfiltered_threshold = (
            self.quality_filtered_reads
            > config.CRITERIA.MIN_FILTERED_READS)
        if qfiltered_threshold:
            if raw_threshold:
                return FLAGS.SUCCESS
            return FLAGS.WARNING
        return FLAGS.DANGER

    @property
    def nanoplot_raw_html_base64(self):
        content_bytes = config.nanoplot_raw_html_path.read_bytes()
        return base64.b64encode(content_bytes).decode()

    @property
    def nanoplot_filtered_html_base64(self):
        content_bytes = config.nanoplot_filtered_html_path.read_bytes()
        return base64.b64encode(content_bytes).decode()

    @property
    def html_file(self):
        return config.run_qc_html_file + '1'


class AbstractResultRows:
    """A result composed of a series of rows with defined columns names."""
    COLUMNS = []

    def __init__(self, rows):
        self.rows = [
            {
                column.lower(): row[column].strip()
                for column in self.COLUMNS
            }
            for row in rows
        ]

    def __len__(self):
        return len(self.rows)

    def __iter__(self):
        return iter(self.rows)

    def __getitem__(self, index):
        return self.rows[index]

    def to_json(self):
        return self.rows


class BlastHits(AbstractResultRows):
    COLUMNS = [
        'Sample_name',
        'qseqid',
        'consensus_sequence',
        'length',
        'reference_title',
        'reference_accession',
        'reference_length',
        'pc_ident',
        'evalue',
        'bitscore',
        'query_coverage',
        'orientation',
        'species',
        'kingdom',
        'full_lineage',
        'read_count',
        'pc_read',
        'pc_depth_30X',
        'target_organism_match',
        'status',
    ]


class BlastHitsPolished(AbstractResultRows):
    COLUMNS = [
        'Sample_name',
        '?',
    ]
