"""Entrypoint for rendering a workflow report."""

import base64
import json
import logging
import os
import csv
from datetime import datetime
from jinja2 import Environment, FileSystemLoader
from pathlib import Path

from . import config


logger = logging.getLogger(__name__)
config = config.Config()

TEMPLATE_DIR = Path(__file__).parent / 'templates'
STATIC_DIR = Path(__file__).parent / 'static'


def render(result_dir: Path):
    """Render to HTML report to the configured output directory."""
    j2 = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    template = j2.get_template('index.html')
    config.read_results(result_dir)
    context = _get_report_context()

    # ! TODO: Remove this
    path = config.output_dir / 'example_report_context.json'
    with path.open('w') as f:
        print(f"Writing report context to {path}")
        json.dump(context, f, indent=2)
    # ! ~~~

    static_files = _get_static_file_contents()
    rendered_html = template.render(**context, **static_files)

    report_path = config.get_report_path()
    with open(report_path, 'w') as f:
        f.write(rendered_html)
    logger.info(f"HTML document written to {report_path}")


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


def _get_report_context() -> dict:
    """Build the context for the report template."""
    return {
        'title': config.REPORT.TITLE,
        'subtitle_html': config.REPORT.SUBTITLE,
        'sample_id': config.sample_id,
        'facility': "Hogwarts",  # ! TODO
        'analyst_name': "John Doe",  # ! TODO
        'start_time': config.start_time.strftime("%Y-%m-%d %H:%M:%S"),
        'end_time': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'wall_time': _get_walltime(),
        'metadata': _get_metadata(),
        'parameters': _get_parameters(),
        'runs_stats': _get_runs_stats(),
    }


def _get_walltime():
    """Return wall time since start of the workflow.
    Returns a dict of hours, minutes, seconds.
    """
    seconds = (datetime.now() - config.start_time).total_seconds()
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return {
        'hours': int(hours),
        'minutes': int(minutes),
        'seconds': int(seconds),
    }


def _get_metadata():
    """Return the metadata as a dict."""
    with config.metadata_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['Sampleid'] == config.sample_id:
                return row


def _get_parameters() -> dict[str, dict[str, str]]:
    """Return the parameters as a dict."""
    with config.parameters_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['Sampleid'] == config.sample_id:
                return row
    return {}


def _get_runs_stats() -> dict:
    """Return the runs stats as a dict.

    Columns:
    - Sample
    - raw_reads
    - quality_filtered_reads
    - percent_quality_filtered
    - raw_reads_flag
    - qfiltered_flag
    """
    with config.runs_stats_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['Sample'] == config.sample_id:
                return row
    return {}
