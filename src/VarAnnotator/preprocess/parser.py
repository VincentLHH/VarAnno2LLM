import argparse
from pathlib import Path
from VarAnnotator.version import __version__

def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="VarAnnotator - Annotation Tool based on VEP",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example：
  varannot input.vcf -o output.vep.vcf \\
    --vep-cache /path/to/vep_cache \\
    --ref-seq /path/to/hg19.fa \\
    --genome GRCh37 \\
    --fork 8 \\
    --verbose

  # remove annotation sources to speed up
  varannot input.vcf -o fast.vcf \\
    --vep-cache ... --ref-seq ... \\
    --no-plugins --no-customs

Menu structure：
  vep_cache/
    ├── vep/             # VEP cache (dir_cache) and plugins (dir_plugins)
    ├── dbNSFP/          # dbNSFP{version}_{assemble}.gz
    ├── spliceAI/        # spliceai_scores.raw.snv.{assemble}.vcf.gz
    ├── dbscSNV/         # dbscSNV{version}_{assemble}.txt.gz
    ├── condel/          # condel config
    ├── clinvar/         # clinvar.vcf.gz
    └── intervar/        # intervar.vcf.gz
"""
    )

    subparsers = parser.add_subparsers(
        dest="command",
        help="Available sub-commands",
        required=True
    )

    parser_pre = subparsers.add_parser(
        "pre",
        help="Preprocess VCF: QC & Normalization"
    )

    parser_pre.add_argument(
        "--verbose", action="store_true", help="output to log"
    )

    parser_pre.add_argument(
        "input",
        type=Path,
        help="Input VCF/VCF.gz (default), or cpra (need --cpra)"
    )

    parser_pre.add_argument(
        "--cpra",
        action="store_true",
        help="Input cpra"
    )

    parser_pre.add_argument(
        "-o", "--output",
        type=Path,
        required=True,
        help="Output VCF.gz"
    )

    parser_pre.add_argument(
        "--ref-seq",
        type=Path,
        required=True,
        help="reference FASTA path"
    )

    parser_pre.add_argument(
        "--min-qual",
        type=float,
        default=30.0,
        help="min QUAL for QC (Default: 30.0)"
    )

    parser_pre.add_argument(
        "--no_qc",
        action="store_true",
        help="pick if not QC (QUAL, FILTER)"
    )

    parser_anno = subparsers.add_parser(
        "anno",
        help="Variant annotation based on VEP"
    )

    parser_anno.add_argument("input", type=Path, help="input VCF path")

    parser_anno.add_argument(
        "-o", "--output", type=Path, required=True, help="output annotated VCF path"
    )

    parser_anno.add_argument(
        "--verbose", action="store_true", help="output to log"
    )

    parser_anno.add_argument(
        "--vep-cache",
        type=Path,
        required=True,
        help="VEP cache path",
    )

    parser_anno.add_argument(
        "--ref-seq",
        type=Path,
        required=True,
        help="reference FASTA path",
    )

    # Optional Params
    parser_anno.add_argument(
        "--genome",
        default="GRCh37",
        choices=["GRCh37", "GRCh38"],
        help="Currently only support: GRCh37",
    )

    parser_anno.add_argument(
        "--output-format",
        default="vcf",
        choices=["vcf", "json", "tab"],
        help="output format (vcf, json, tab)",
    )

    parser_anno.add_argument(
        "--vep-path",
        default="conda run -n vep vep",
        help="VEP cmd (Default: conda run "
             "-n vep vep)",
    )

    parser_anno.add_argument(
        "--fork",
        type=int,
        default=8,
        help="forks (default:8)",
    )

    parser_anno.add_argument(
        "--no-overwrite",
        action="store_true",
        help="do not overwrite existing annotated VCF",
    )

    # Plugin switch
    parser_anno.add_argument(
        "--no-plugins",
        action="store_true",
        help="ban all VEP plugins (dbNSFP, SpliceAI, dbscSNV, Condel)",
    )

    parser_anno.add_argument(
        "--no-customs",
        action="store_true",
        help="ban all VEP customs (ClinVar, InterVar, OMIM, HGMD)",
    )

    parser_anno.add_argument(
        "--version", action="version", version=f"VarAnnotator {__version__}"
    )

    parser_prompt = subparsers.add_parser(
        "prompt",
        help="Prompting annotation data for LLM"
    )

    parser_prompt.add_argument(
        "input", type=Path, help="Input Annotated File path (Tab separated)"
    )

    parser_prompt.add_argument(
        "--gene-list",
        type=Path,
        required=True,
        help="gene list (for RAG) path",
    )

    parser_prompt.add_argument(
        "--white-list",
        type=Path,
        default=None,
        help="white list variants",
    )

    parser_prompt.add_argument(
        "--verbose", action="store_true", help="output to log"
    )

    parser_prompt.add_argument(
        "--pheno-file",
        type=Path,
        default=None,
        help="pheno list path",
    )

    parser_prompt.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="output jsonl folder (*.jsonl)",
    )

    parser_prompt.add_argument(
        "--af-filter",
        action="store_true",
        help="use af filter (<0.01)",
    )

    parser_prompt.add_argument(
        "--no-whitelist",
        action="store_true",
        help="not use white gene list",
    )

    parser_prompt.add_argument(
        "--label",
        choices=["test", "positive", "negative"],
        default="test",
        help="prompt label",
    )

    parser_prompt.add_argument(
        "--consequence-filter",
        action="store_true",
        help="keep protein-altering and splice-related consequences only",
    )

    parser_prompt.add_argument(
        "--omim-filter",
        action="store_true",
        help="keep variants with a non-empty OMIM disease annotation",
    )

    parser_prompt.add_argument(
        "--max-variant-length",
        type=int,
        default=None,
        help="exclude alleles whose REF or ALT length is at least this value",
    )


    return parser
