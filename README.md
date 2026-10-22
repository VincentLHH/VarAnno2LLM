# VarAnno2LLM

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

**VarAnno2LLM** is an end-to-end variant annotation and prompt generation pipeline designed to transform raw genomic variants (VCF or CPRA format) into LLM-friendly JSONL files for training or inference. The toolkit integrates three core functionalities:

1. **Preprocessing**: QC, normalization, and aggregation of sparse CPRA/phenotype data.
2. **Annotation**: Comprehensive functional annotation powered by [Ensembl VEP](https://grch37.ensembl.org/info/docs/tools/vep/index.html) with support for plugins and custom databases.
3. **Prompt Engineering**: Conversion of annotated variants into structured JSONL records embedded with customizable prompt templates for LLMs.

> Contact: liuhaohan@genomics.cn

---

## Core Features

### 1. Preprocessing (`pre`)
- Accepts VCF/VCF.gz or CPRA input
- Performs quality control (QUAL/FILTER filtering), left-alignment/normalization and multiallelic splitting (via `bcftools norm`)
- Option to skip QC for pre-cleaned data

### 2. Annotation (`anno`)
- Runs VEP with rich plugin support:
  - **dbNSFP**: pathogenicity & conservation scores
  - **SpliceAI**: splice-altering predictions
  - **dbscSNV**: splice site variants
  - **Condel**: consensus deleteriousness
- Supports custom annotation tracks:
  - ClinVar, InterVar, OMIM, HGMD (via VEP `--custom`)
- Flexible toggles: disable plugins or custom annotations for faster runs

### 3. Prompt Generation (`prompt`)
- Converts VEP-annotated TSV into **LLM-ready JSONL**
- Supports:
  - Gene-list membership annotation for RAG
  - A known-variant whitelist for allele-frequency filter exemptions
  - Phenotype integration (HPO or custom phenotypes)
  - Allele frequency filtering (AF < 0.01)
  - Optional protein-altering/splice consequence filtering
  - Optional OMIM-annotated-gene filtering
  - Optional allele-length filtering
  - Label assignment (`test`, `positive`, `negative`)
- Embeds configurable prompt templates (editable in source)

---

## Cache Directory Structure

```bash
vep_cache/
├── vep/             # VEP cache + plugins
├── dbNSFP/          # dbNSFP{version}_{assembly}.gz
├── spliceAI/        # spliceai_scores.raw.snv.{assembly}.vcf.gz
├── dbscSNV/         # dbscSNV{version}_{assembly}.txt.gz
├── condel/          # Condel configuration files
├── clinvar/         # clinvar.vcf.gz
└── intervar/        # intervar.vcf.gz
```

---

## Quick Start

You can run **VarAnno2LLM** locally in minutes:

```bash
# 1. Clone the repository
git clone https://github.com/VincentLHH/VarAnno2LLM.git
cd VarAnno2LLM

# 2. Install in editable mode (requires Python ≥3.10)
pip install -e .

# 3. Run any subcommand
varannot --help
```

Run the test suite after the editable install:

```bash
python -m pytest -q
```

If the package has not yet been installed into the active environment, use:

```bash
PYTHONPATH=src python -m pytest -q
```

VIPER-style prompt generation uses the following candidate filters:

```bash
varannot prompt annotated.tsv \
  --gene-list gene_list.xlsx \
  -o prompts \
  --af-filter \
  --consequence-filter \
  --omim-filter \
  --max-variant-length 1000
```

## External resources and versioning

VarAnno2LLM calls external VEP and bcftools executables. For reproducible
results, record the VEP release and cache version, GRCh37 FASTA source,
dbNSFP, SpliceAI, dbscSNV, ClinVar, InterVar, HGMD and OMIM snapshot dates in
the run metadata. Archive the environment or container together with the exact
source commit used for a publication.


