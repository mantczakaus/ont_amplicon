import os
from jinja2 import Environment, FileSystemLoader
from pathlib import Path

# Paths
PWD = Path(__file__).parent
TEMPLATE_DIR = PWD / 'templates'
TEMPLATE_NAME = "bam-viewer.html"
OUTPUT_DIR = PWD / 'output'
DATA_DIR = Path('/home/cameron/Downloads/data/bam')
SAMPLE_ID = 'BC11_MOL-010c_LMJ1576_C'
BAM_FILE_PATH = DATA_DIR / "BC11_MOL-010c_LMJ1576_C_aln.sorted.bam"
BAI_FILE_PATH = DATA_DIR / "BC11_MOL-010c_LMJ1576_C_aln.sorted.bam.bai"
FASTA_FILE_PATH = (
    DATA_DIR / "BC11_MOL-010c_LMJ1576_C_final_polished_consensus_match.fasta")
FAI_FILE_PATH = (
    DATA_DIR /
    "BC11_MOL-010c_LMJ1576_C_final_polished_consensus_match.fasta.fai")
OUTPUT_HTML_PATH = OUTPUT_DIR / TEMPLATE_NAME


def file_to_js_array(path):
    byte_array = path.read_bytes()
    return "[" + ",".join(str(b) for b in byte_array) + "]"


# Check if BAM and BAI files exist
if not os.path.exists(BAM_FILE_PATH) or not os.path.exists(BAI_FILE_PATH):
    raise FileNotFoundError("Both BAM and BAI files are required.")

j2 = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
template = j2.get_template(TEMPLATE_NAME)
context = {
    k: file_to_js_array(v)
    for k, v in [
        ('bam_binary_arr', BAM_FILE_PATH),
        ('bai_binary_arr', BAI_FILE_PATH),
        ('fasta_binary_arr', FASTA_FILE_PATH),
        ('fai_binary_arr', FAI_FILE_PATH),
    ]
}
rendered_html = template.render(sample_id=SAMPLE_ID, **context)

# Write output HTML
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
with open(OUTPUT_HTML_PATH, "w") as f:
    f.write(rendered_html)

print(f"BAM Viewer HTML generated: {OUTPUT_HTML_PATH}")
