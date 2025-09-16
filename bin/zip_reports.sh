#!/usr/bin/bash

# Zip all workflow reports, bam viewers, QC reports for easy download

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --dir) DIR="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$DIR" ]]; then
    echo "Usage: $0 --dir <directory>"
    exit 1
fi

cd "$DIR"
mkdir reports

for d in */; do
    sample=$(basename "$d")
    if [[ -d "$sample/07_html_report" ]]; then
        mkdir -p "reports/$sample"
        cp "$sample"/07_html_report/*_report.html "reports/$sample"
        cp "$sample"/07_html_report/run_qc_report.html "reports/$sample"
        cp "$sample"/07_html_report/*_bam-alignment.html "reports/$sample"
        cp "$sample"/03_polishing/*_final_polished_consensus.fasta "reports/$sample"
    fi
done

zip -r reports.zip reports/ > /dev/null

echo "Primary workflow outputs have been zipped"
echo 'You can download all workflow reports and sequences by clicking the following files under "Results" tab:'
echo "- reports.zip"
echo "- sequences.zip"
