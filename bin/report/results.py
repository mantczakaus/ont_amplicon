"""Define specific results from the analysis."""

import base64
from Bio import SeqIO

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
        return config.run_qc_html_file


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
    COLUMN_LABELS = [
        ('sample_name', 'Sample ID'),
        ('qseqid', 'Query ID'),
        ('consensus_seq', None),
        ('sgi', None),
        ('sacc', 'Hit accession'),
        ('alignment_length', 'Alignment length'),
        ('nident', None),
        ('pident', 'Identity %'),
        ('mismatch', None),
        ('gaps', 'Gaps'),
        ('gapopen', None),
        ('qstart', 'Query start'),
        ('qend', 'Query end'),
        ('qlen', None),
        ('sstart', 'Subject start'),
        ('send', 'Subject end'),
        ('slen', 'Subject length'),
        ('sstrand', None),
        ('evalue', 'E-value'),
        ('bitscore', 'Bitscore'),
        ('qcovhsp', None),
        ('stitle', 'Subject title'),
        ('staxids', None),
        ('qseq', None),
        ('sseq', None),
        ('sseqid', 'Subject ID'),
        ('qcovs', None),
        ('qframe', None),
        ('sframe', None),
        ('species', None),
        ('sskingdoms', 'Kingdom'),
        ('FullLineage', None),
        ('target_organism_match', 'Target organism match'),
        ('n_read_cont_cluster', None),
        ('query_match_length', None),
        ('qseq_mapping_read_count', 'Reads mapped'),
        ('qseq_mean_depth', 'Read depth'),
        ('qseq_pc_mapping_read', 'Reads mapped (%)'),
        ('qseq_pc_depth_30X', '30X depth coverage (%)'),
        ('30X_DEPTH_FLAG', None),
        ('MAPPED_READ_COUNT_FLAG', None),
        ('TARGET_ORGANISM_FLAG', None),
        ('TARGET_SIZE_FLAG', None),
    ]
    COLUMNS = [c for c, _ in COLUMN_LABELS]
    COLUMNS_TO_DISPLAY = [c for c in COLUMN_LABELS if c[1] is not None]


class BlastHitsPolished(AbstractResultRows):
    COLUMNS = [
        'qseqid',
        'sgi',
        'sacc',
        'length',
        'nident',
        'pident',
        'mismatch',
        'gaps',
        'gapopen',
        'qstart',
        'qlen',
        'qend',
        'sstart',
        'send',
        'slen',
        'sstrand',
        'evalue',
        'bitscore',
        'qcovhsp',
        'stitle',
        'staxids',
        'qseq',
        'sseq',
        'sseqid',
        'qcovs',
        'qframe',
        'sframe',
        'species',
        'sskingdoms',
    ]


class ConsensusFASTA:

    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        self.records = list(SeqIO.parse(fasta_path, 'fasta'))

    def __len__(self):
        return len(self.records)

    def __iter__(self):
        return iter(self.records)

    def to_json(self):
        return {
            seq.id: str(seq.seq)
            for seq in self.records
        }
