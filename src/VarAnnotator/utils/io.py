import pandas as pd
from pathlib import Path
from typing import Iterator, Dict, Any
from VarAnnotator.utils.logging import setup_logger
import logging
setup_logger()
logger = logging.getLogger(__name__)

def get_anno_dict(filepath: str) -> Iterator[Dict[str, Any]]:
    """
    Parse VEP-style TSV file from a file path.

    Args:
        filepath (str): Path to the VEP TSV file.

    Yields:
        dict: Annotation record for each variant.
    """
    header = None

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("##"):
                continue

            if line.startswith("#"):
                header = line[1:].strip().split("\t")
                continue

            if header is None:
                logger.error("Header not found. The file must have a # header line.")

            fields = line.split("\t")
            if len(fields) != len(header):
                logger.error(
                    f"Field count mismatch: {len(fields)} fields, {len(header)} headers.\n"
                    f"Line: {line}"
                )

            record = {}
            for key, value in zip(header, fields):
                value = value.strip()
                if not value or value == "-":
                    record[key] = None
                elif "," in value and key != "OMIM" and key != 'Pheno_PHENOTYPE':
                    values = [v.strip() for v in value.split(",")]
                    cleaned = [v for v in values if v not in (".", "")]
                    record[key] = cleaned if cleaned else None
                else:
                    record[key] = value if value != "." else None

            yield record

def get_gene_set(gene_list_file: Path) -> set:
    df = pd.read_csv(gene_list_file)
    gene_set = set(df['Gene Symbol'])
    return gene_set

def get_white_list(white_list_file: Path) -> set:
    df = pd.read_csv(white_list_file)
    white_set = set(df['cpra'])
    return white_set

def get_pheno_info(pheno_file: Path) -> dict:
    df = pd.read_csv(pheno_file)
    return pd.Series(df.phenotype.values, index=df.sample_name).fillna("").astype(str).to_dict()
