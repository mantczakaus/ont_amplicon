from jinja2 import Environment, FileSystemLoader
from pathlib import Path

TEMPLATE_DIR = Path(__file__).parent / 'templates'
TEMPLATE_NAME = "bam-viewer.html"


def file_to_js_array(path):
    byte_array = path.read_bytes()
    return "[" + ",".join(str(b) for b in byte_array) + "]"


def render_html(src_dir, output_dir):
    bam_file = list(src_dir.glob("*.bam"))[0]
    bai_file = list(src_dir.glob("*.bai"))[0]
    fasta_file = list(src_dir.glob("*consensus_match.fasta"))[0]
    sample_id = bam_file.stem.split('_aln.sorted')[0]

    j2 = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    template = j2.get_template(TEMPLATE_NAME)
    context = {
        k: file_to_js_array(v)
        for k, v in [
            ('bam_binary_arr', bam_file),
            ('bai_binary_arr', bai_file),
            ('fasta_binary_arr', fasta_file),
        ]
    }
    rendered_html = template.render(sample_id=sample_id, **context)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "bam.html"
    path.write_text(rendered_html)
    print(f"BAM Viewer HTML generated: {path}")
