"""Render a workflow report from output files."""

import base64
import json
import logging
import os
import csv
from datetime import datetime
from jinja2 import Environment, FileSystemLoader
from pathlib import Path

from . import config
from .results import (
    BlastHits,
    ConsensusFASTA,
    Metadata,
    RunQC,
)
from .utils import serialize
from .bam import render_bam_html

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
config = config.Config()

TEMPLATE_DIR = Path(__file__).parent / 'templates'
STATIC_DIR = Path(__file__).parent / 'static'


def render(result_dir: Path, samplesheet: Path):
    """Render to HTML report to the configured output directory."""
    config.load(result_dir)
    render_bam_html()
    j2 = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    template = j2.get_template('index.html')
    context = _get_report_context(samplesheet)

    # ! TODO: Remove this
    path = config.result_dir / 'example_report_context.json'
    with path.open('w') as f:
        logger.info(f"Writing report context to {path}")
        json.dump(context, f, indent=2, default=serialize)
    # ! ~~~

    static_files = _get_static_file_contents()
    rendered_html = template.render(**context, **static_files)

    with open(config.report_path, 'w') as f:
        f.write(rendered_html)
    logger.info(f"HTML document written to {config.report_path}")


def _get_static_file_contents():
    """Return the static files content as strings."""
    static_files = {}
    for root, _, files in os.walk(STATIC_DIR):
        root = Path(root)
        if root.name == 'css':
            static_files['css'] = [
                f'/* {f} */\n' + (root / f).read_text()
                for f in files
            ]
        elif root.name == 'js':
            static_files['js'] = [
                f'/* {f} */\n' + (root / f).read_text()
                for f in files
            ]
        elif root.name == 'img':
            static_files['img'] = {
                f: _get_img_src(root / f)
                for f in files
            }
    return {'static': static_files}


def _get_img_src(path):
    """Return the base64 encoded image source as an HTML img src property."""
    ext = path.suffix[1:]
    return (
        f"data:image/{ext};base64,"
        + base64.b64encode(path.read_bytes()).decode()
    )


def _get_report_context(samplesheet) -> dict:
    """Build the context for the report template."""
    blast_hits = _get_blast_hits()
    consensus_fasta = ConsensusFASTA(config.consensus_fasta_path)
    consensus_match_fasta = ConsensusFASTA(config.consensus_match_fasta_path)
    return {
        'title': config.REPORT.TITLE,
        'subtitle_html': config.REPORT.SUBTITLE,
        'sample_id': config.sample_id,
        'facility': "Hogwarts",  # ! TODO
        'analyst_name': "John Doe",  # ! TODO
        'start_time': _get_start_time(),
        'end_time': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'wall_time': _get_walltime(),
        'metadata': _get_metadata(samplesheet),
        'parameters': {},  # _get_parameters(),  # TODO: doesn't exist yet
        'run_qc': _get_run_qc(),
        'bam_html_file': config.bam_html_output_path.name,
        'consensus_blast_hits': blast_hits,
        'consensus_blast_stats': {
            'percent': round(100 * len(blast_hits) / len(consensus_fasta)),
            'count': len(blast_hits),
        },
        'consensus_fasta': consensus_fasta,
        'consensus_match_fasta': consensus_match_fasta,
    }


def _get_start_time():
    if not config.start_time:
        return None
    return config.start_time.strftime("%Y-%m-%d %H:%M:%S")


def _get_walltime():
    """Return wall time since start of the workflow.
    Returns a dict of hours, minutes, seconds.
    """
    if not config.start_time:
        return None
    seconds = (datetime.now() - config.start_time).total_seconds()
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return {
        'hours': str(int(hours)).zfill(2),
        'minutes': str(int(minutes)).zfill(2),
        'seconds': str(int(seconds)).zfill(2),
    }


def _get_metadata(samplesheet):
    """Return the metadata as a dict."""
    with open(samplesheet) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['sampleid'] == config.sample_id:
                return Metadata(row)


def _get_parameters() -> dict[str, dict[str, str]]:
    """Return the parameters as a dict."""
    with config.parameters_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['Sampleid'] == config.sample_id:
                return row
    return {}


def _get_run_qc() -> dict:
    """Return the runs stats as a dict.

    Columns:
    - Sample
    - raw_reads
    - quality_filtered_reads
    - percent_quality_filtered
    - raw_reads_flag
    - qfiltered_flag
    """
    with config.run_qc_path.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['Sample'] == config.sample_id:
                return RunQC(row)
    return {}


def _get_blast_hits() -> BlastHits:
    """Return the blast hits as a dict."""
    with config.blast_hits_path.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        return BlastHits(reader)
