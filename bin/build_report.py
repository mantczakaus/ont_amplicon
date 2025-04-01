#!/usr/bin/env python

"""Build the HTML report from output data."""

import argparse

from report import report
from report.utils import existing_path


def main():
    """Parse the command line arguments and build the report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", type=str)
    parser.add_argument(
        '--result_dir',
        type=existing_path,
        help="The directory containing the output data.",
    )
    args = parser.parse_args()
    report.render(args.result_dir, args.samplesheet)


if __name__ == '__main__':
    main()
