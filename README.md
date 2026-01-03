# VarAnno2LLM

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**VarAnno2LLM** is an end-to-end variant annotation and prompt generation pipeline designed to transform raw genomic variants (VCF or CPRA format) into LLM-friendly JSONL files for training or inference. The toolkit integrates three core functionalities:

1. **Preprocessing**: QC, normalization, and aggregation of sparse CPRA/phenotype data with genotype-based filtering (e.g., remove high-frequency variants, apply whitelist).
2. **Annotation**: Comprehensive functional annotation powered by [Ensembl VEP](https://grch37.ensembl.org/info/docs/tools/vep/index.html) with support for plugins and custom databases.
3. **Prompt Engineering**: Conversion of annotated variants into structured JSONL records embedded with customizable prompt templates for LLMs.

> Contact: [liuhaohan@genomics.cn](mailto:liuhaohan@genomics.cn)

---

## Core Features

### 1. Preprocessing (`pre`)
- Accepts VCF/VCF.gz or CPRA input
- Performs quality control (QUAL/FILTER filtering) and left-alignment/normalization (via `bcftools norm`)
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
  - Gene whitelist filtering (for RAG or targeted analysis)
  - Phenotype integration (HPO or custom phenotypes)
  - Allele frequency filtering (AF < 0.01)
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

