"""Entrypoint for rendering a workflow report."""

import base64
import csv
import json
import logging
import os
from datetime import datetime
from jinja2 import Environment, FileSystemLoader
from pathlib import Path

from .outcomes import DetectedTaxon
from ..utils import config
from ..utils.flags import Flag, FLAGS  # , TARGETS

logger = logging.getLogger(__name__)
config = config.Config()

TEMPLATE_DIR = Path(__file__).parent / 'templates'
STATIC_DIR = Path(__file__).parent / 'static'


def render(query_ix):
    """Render to HTML report to the configured output directory."""
    j2 = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    template = j2.get_template('index.html')
    context = _get_report_context(query_ix)

    # ! TODO: Remove this
    path = config.output_dir / 'example_report_context.json'
    with path.open('w') as f:
        from src.utils import serialize
        print(f"Writing report context to {path}")
        json.dump(context, f, default=serialize, indent=2)
    # ! ~~~

    static_files = _get_static_file_contents()
    rendered_html = template.render(**context, **static_files)

    report_path = config.get_report_path(query_ix)
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


def _get_report_context(query_ix):
    """Build the context for the report template."""
    query_fasta_str = config.read_query_fasta(query_ix).format('fasta')
    return {
        'title': config.REPORT.TITLE,
        'facility': "Hogwarts",  # ! TODO
        'analyst_name': "John Doe",  # ! TODO
        'start_time': config.start_time.strftime("%Y-%m-%d %H:%M:%S"),
        'end_time': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'wall_time': _get_walltime(),
        'metadata': _get_metadata(query_ix),
        'config': config,
        'input_fasta': query_fasta_str,
        'conclusions': _draw_conclusions(query_ix),
        'hits': config.read_blast_hits_json(query_ix)['hits'],
        'candidates': _get_candidates(query_ix),
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


def _get_metadata(query_ix):
    """Return mock metadata for the report."""
    return [  # TODO: parse from metadata.csv
        {
            'name': 'Sample ID',
            'value': 'LC438549.1',
        },
        {
            'name': 'Locus',
            'value': 'COI',
        },
        {
            'name': 'Preliminary morphology ID',
            'value': 'Aphididae',
        },
        {
            'name': 'Taxa of interest',
            'value': [
                'Myzus persicae',
            ],
        },
        {
            'name': 'Country of origin',
            'value': 'Ecuador',
        },
        {
            'name': 'Host/commodity of origin',
            'value': 'Cut flower Rosa',
        },
        {
            'name': 'Comments',
            'value': 'Lorem ipsum dolor sit amet',
        },
    ]


def _draw_conclusions(query_ix):
    """Determine conclusions from outputs flags and files."""
    flags = Flag.read(query_ix)
    return {
        'flags': {
            ix: flag.to_json()
            for ix, flag in flags.items()
        },
        'summary': {
            'result': _get_taxonomic_result(query_ix, flags),
            'pmi': _get_pmi_result(flags),
            'toi': _get_toi_result(query_ix, flags),
        }
    }


def _get_taxonomic_result(query_ix, flags):
    """Determine the taxonomic result from the flags."""
    path = config.get_query_dir(query_ix) / config.TAXONOMY_ID_CSV
    flag_1 = flags[FLAGS.POSITIVE_ID]
    explanation = (f"Flag {FLAGS.POSITIVE_ID}{flag_1.value}: "
                   + flag_1.explanation())
    if flag_1.value == FLAGS.A:
        with path.open() as f:
            reader = csv.DictReader(f)
            hit = next(reader)
        return {
            'confirmed': True,
            'explanation': explanation,
            'species': hit['species'],
        }
    return {
        'confirmed': False,
        'explanation': explanation,
        'species': None,
    }


def _get_pmi_result(flags):
    """Determine the preliminary ID confirmation from the flags."""
    flag_1 = flags[FLAGS.POSITIVE_ID]
    if flag_1.value != FLAGS.A:
        return {
            'confirmed': False,
            'explanation': "Inconclusive taxonomic identity (Flag"
                           f" {FLAGS.POSITIVE_ID}{flag_1.value})",
            'bs-class': None,
        }
    flag_7 = flags[FLAGS.PMI]
    if flag_7.value == FLAGS.A:
        return {
            'confirmed': True,
            'explanation': flag_7.explanation(),
            'bs-class': 'success',
        }
    return {
        'confirmed': False,
        'explanation': flag_7.explanation(),
        'bs-class': 'danger',
    }


def _get_toi_result(query_ix, flags):
    """Determine the taxa of interest detection from the flags."""
    query_dir = config.get_query_dir(query_ix)
    path = query_dir / config.TOI_DETECTED_CSV
    with path.open() as f:
        reader = csv.DictReader(f)
        detected_tois = [
            DetectedTaxon(*[
                row.get(colname)
                for colname in config.OUTPUTS.TOI_DETECTED_HEADER
            ])
            for row in reader
            if row.get(config.OUTPUTS.TOI_DETECTED_HEADER[1])
        ]
    flag_2 = flags[FLAGS.TOI]
    criteria_2 = f"Flag {flag_2}: {flag_2.explanation()}"

    # TODO
    # flag_5_1 = flags[FLAGS.DB_COVERAGE_TARGET]
    # flag_5_2 = flags[FLAGS.DB_COVERAGE_RELATED]
    # flag_5_1_value = flag_5_1.value_for_target(TARGETS.TOI, max_only=True)
    # criteria_5_1 = (
    #     f"Flag {flag_5_1} for taxa of concern:"
    #     f" {flag_5_1.explanation(flag_5_1_value)}")
    # flag_5_2_value = flag_5_2.value_for_target(TARGETS.TOI, max_only=True)
    # criteria_5_2 = (
    #     f"Flag {flag_5_2} for taxa of concern:"
    #     f" {flag_5_2.explanation(flag_5_2_value)}")
    return {
        'detected': detected_tois,
        'criteria': [
            {
                'message': criteria_2,
                'level': flag_2.get_level(),
                'bs-class': flag_2.get_bs_class(),
            },
            # { # TODO
            #     'message': criteria_5_1,
            #     'level': flag_5_1.get_level(flag_5_1_value),
            #     'bs-class': flag_5_1.get_bs_class(flag_5_1_value),
            # },
            # {
            #     'message': criteria_5_2,
            #     'level': flag_5_2.get_level(flag_5_2_value),
            #     'bs-class': flag_5_2.get_bs_class(flag_5_2_value),
            # },
        ],
        'ruled_out': (
            flag_2.value == FLAGS.A
            # and flag_5_1_value == FLAGS.A  # TODO
            # and flag_5_2_value == FLAGS.A  # TODO
        ),
        'bs-class': 'success' if flag_2.value == FLAGS.A else 'danger',
    }


def _get_candidates(query_ix):
    """Read data for the candidate hits/taxa."""
    flags = Flag.read(query_ix)
    query_dir = config.get_query_dir(query_ix)
    with open(query_dir / config.CANDIDATES_JSON) as f:
        candidates = json.load(f)
    candidates['fasta'] = {
        seq.id: seq.format("fasta")
        for seq in config.read_fasta(query_dir / config.CANDIDATES_FASTA)
    }
    candidates['strict'] = (
        flags[FLAGS.POSITIVE_ID].value
        not in (FLAGS.D, FLAGS.E))
    return candidates


if __name__ == '__main__':
    query_ix = 0
    render(query_ix)
