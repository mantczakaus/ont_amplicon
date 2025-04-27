"""Define specific results from the analysis."""

import base64
import csv
from Bio import SeqIO

from .config import Config

config = Config()


def _csv_to_dict(csv_path, index_col='colname'):
    with open(csv_path) as f:
        reader = csv.reader(f)
        next(reader)  # skip header line
        ordered_colnames = [
            row[0].strip()
            for row in reader
        ]
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        data = {
            row[index_col].strip(): {
                colname: value
                for colname, value in row.items()
            }
            for row in reader
        }
        return {
            colname: data[colname]
            for colname in ordered_colnames
        }


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
    COLUMN_METADATA = {}

    def __init__(self, rows):
        self.rows = self._parse_rows(rows)

    def __len__(self):
        return len(self.rows)

    def __iter__(self):
        return iter(self.rows)

    def __getitem__(self, index):
        return self.rows[index]

    def _parse_rows(self, rows):
        return [
            {
                colname: (
                    self._cast(
                        row[colname].strip(),
                        self.COLUMN_METADATA.get(colname, {}).get('type'),
                    )
                )
                for colname in self.COLUMNS
            }
            for row in rows
        ]

    def _cast(self, value, type_str):
        if not type_str:
            return value
        if type_str == 'int':
            num = int(float(value))
            return f"{num:,}"
        if type_str == 'float':
            return float(value)
        if type_str == 'scientific' and 'e' in value:
            num = float(value)
            return f"{num:.2e}"
        return value

    def to_json(self):
        return self.rows


class BlastHits(AbstractResultRows):
    COLUMN_METADATA = _csv_to_dict(config.SCHEMA.BLAST_HITS_FIELD_CSV)
    COLUMNS = list(COLUMN_METADATA.keys())

    def __init__(self, *args):
        super().__init__(*args)
        self.columns_display = [
            c for c in self.COLUMN_METADATA
            if self.COLUMN_METADATA[c]['label']
        ]
        self.columns_primary_display = [
            c for c in self.columns_display
            if self.COLUMN_METADATA[c]['primary_display']
        ]
        self.rows = self.set_bs_class()

    def set_bs_class(self):
        def _get_bs_class(row):
            flags = (
                row['30X_DEPTH_FLAG'],
                row['MAPPED_READ_COUNT_FLAG'],
            )
            if 'ORANGE' in flags:
                return 'warning'
            if 'RED' in flags:
                return 'danger'
            if 'GREY' in flags:
                return 'secondary'
            return 'success'

        return [
            {
                **row,
                'bs_class': _get_bs_class(row),
            }
            for row in self.rows
        ]


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
        'broad_taxonomic_category',
    ]


class ConsensusFASTA:

    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        self.records = list(SeqIO.parse(fasta_path, 'fasta'))

    def __len__(self):
        return len(self.records)

    def __iter__(self):
        return iter(self.records)

    def __str__(self):
        return '\n\n'.join(
            f">{seq.id}\n"
            + self._wrap(seq.seq)
            for seq in self.records
        )

    def _wrap(self, seq, width=80):
        seq_str = str(seq)
        return '\n'.join(
            seq_str[i:i + width]
            for i in range(0, len(seq), width)
        )

    def to_json(self):
        return {
            seq.id: str(seq.seq)
            for seq in self.records
        }
