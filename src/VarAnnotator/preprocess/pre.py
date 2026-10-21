from pathlib import Path
from typing import Union, Dict, Optional
import subprocess
from VarAnnotator.utils.logging import setup_logger
import logging

setup_logger()
logger = logging.getLogger(__name__)

def cpra_to_vcf(
        cpra_input: Union[str, Path],
        output_vcf: Path,
        sample_name: str = "SAMPLE1"
) -> bool:
    """
    Construct a VCF from CPRA input (CPRA, phenotype, zygosity as tab-separated fields).

    Zygosity can be: Het, Hom, Hemi, or empty (treated as Het with a warning).
    QUAL is set to 100 for all variants.
    GT is determined by zygosity:
        - Hom -> 1/1
        - Het -> 0/1
        - Hemi -> 0/1
        - empty -> 1/1 (default to Hom)

    The output VCF is bgzipped and indexed with tabix.

    Args:
        cpra_input: Path to input file or a string of one cpra entry.
        output_vcf: Output VCF path (will be .vcf.gz).
        sample_name: Sample name in VCF.

    Returns:
        bool: True if VCF generated successfully.
    """
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    if isinstance(cpra_input, (str, Path)) and Path(cpra_input).exists():
        with open(cpra_input, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
    else:
        # Assume it's a single string input
        lines = [str(cpra_input)] if cpra_input else []

    # Map zygosity to GT (genotype)
    zygosity_to_gt: Dict[str, str] = {
        'Hom': '1/1',
        'Het': '0/1',
        'Hemi': '0/1',
        '': '0/1', # default to Het and warn
        '.': '0/1'
    }
    rraw_vcf = output_vcf.with_suffix(".vcf")
    compressed_vcf = str(rraw_vcf) + ".gz"
    chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    try:
        raw_vcf = f"unsorted_{rraw_vcf}"
        with open(raw_vcf, 'w') as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("##source=cpra_to_vcf\n")
            for chrom in chromosomes:
                f.write(f"##contig=<ID={chrom}>\n")
            f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">\n')
            f.write('##INFO=<ID=ORIG_CPRA,Number=1,Type=String,Description="Original CPRA notation (chrom-pos-ref-alt)">\n')
            f.write('##INFO=<ID=PHENOTYPE,Number=1,Type=String,Description="Phenotype info">\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths for the ref and alt alleles">\n')
            f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
            f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")
            ref_depth = 50
            alt_depth = 50
            dp = ref_depth + alt_depth
            ad = f"{ref_depth},{alt_depth}"
            for line_num, line in enumerate(lines, 1):
                parts = line.split('\t')
                if len(parts) < 1:
                    logging.warning(f"Skipping empty line: {line_num}")
                    continue

                # Handle both 1-column (no zygo) and 2-column input
                cpra = parts[0].strip()
                ori_phenotype = parts[1].strip() if len(parts) >= 2 else ''
                phenotype = ori_phenotype.replace(';', '_').replace(',', '_').replace('\t', '' ).replace('|','_').strip()
                zygosity = parts[2].strip() if len(parts) >= 3 else ''
                zygosity = zygosity.capitalize()  # Normalize: het -> Het, etc.

                if zygosity not in zygosity_to_gt:
                    logger.warning(f"Unknown zygosity '{zygosity}' at line {line_num}, treating as Het")
                    zygosity = ''
                if not zygosity:
                    logger.warning(
                        f"Missing zygosity at line {line_num}; using Het for VCF conversion"
                    )
                try:
                    chrom, pos_str, ref, alt = cpra.replace('-', '\t').split('\t')
                    pos = int(pos_str)
                except Exception as e:
                    logger.error(f"Invalid Cpra format at line {line_num}: {cpra}, error: {e}")
                    return False
                gt = zygosity_to_gt[zygosity]
                info_field = f"DP={dp};ORIG_CPRA={chrom}-{pos}-{ref}-{alt};PHENOTYPE={phenotype}"
                format_field = "GT:AD:DP"
                sample_data = f"{gt}:{ad}:{dp}"
                # Fixed QUAL=100, FILTER=PASS
                f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t100\tPASS\t{info_field}\t{format_field}\t{sample_data}\n")

        logger.info(f"Wrote raw VCF: {raw_vcf}")

        # Sort with bcftools
        subprocess.run(["bcftools", "sort", str(raw_vcf), "-o", str(rraw_vcf)], check=True)
        logger.info(f"Sorted VCF: {rraw_vcf}")

        # Compress with bgzip
        subprocess.run(["bgzip", "-f", str(rraw_vcf)], check=True)
        logger.info(f"Compressed VCF: {compressed_vcf}")

        # Index with tabix
        subprocess.run(["tabix", "-f", "-p", "vcf", compressed_vcf], check=True)
        logger.info(f"Indexed VCF: {compressed_vcf}.tbi")

        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"External tool failed: {e}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error during VCF generation: {e}")
        return False


def qc_vcf(
        input_vcf: Path,
        output_vcf: Path,
        ref_dir: Path,
        min_qual: float = 30.0,
        pass_only: bool = True,
        norm: bool = True,
) -> bool:
    """
    Apply QC filter and optional normalization to VCF.

    Filters by QUAL and PASS status, left-aligns and decomposes with bcftools norm.
    Outputs compressed VCF (.vcf.gz) and creates tabix index.

    Args:
        input_vcf: Input VCF file (can be gzipped).
        output_vcf: Output path (without .gz; will generate .vcf.gz and .tbi).
        ref_dir: Reference genome FASTA file.
        min_qual: Minimum QUAL threshold.
        pass_only: Only keep variants with FILTER==PASS or unfiltered (".").
        norm: Whether to run bcftools norm for left-alignment and multi-allelic decomposition.

    Returns:
        bool: True if success.
    """
    input_vcf = Path(input_vcf)
    output_vcf = Path(output_vcf)
    ref_dir = Path(ref_dir)

    if not input_vcf.exists():
        logger.error(f"Input VCF not found: {input_vcf}")
        return False

    if not ref_dir.exists():
        logger.error(f"Reference FASTA not found: {ref_dir}")
        return False

    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    compressed_vcf = str(output_vcf) + ".gz"

    view_cmd = [
        "bcftools", "view",
        "-Ou",
        str(input_vcf)
    ]

    if pass_only:
        filter_expr = f'QUAL >= {min_qual} && (FILTER=="PASS" || FILTER==".")'
        view_cmd += ["-i", filter_expr]

    try:
        result_view = subprocess.run(
            view_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True
        )

        if norm:
            norm_cmd = [
                "bcftools", "norm", "-c", "s", "-m", "-any",
                "-f", str(ref_dir),
                "-Oz",
                "-o", compressed_vcf
            ]
            result_norm = subprocess.run(
                norm_cmd,
                input=result_view.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True
            )
            if result_norm.returncode != 0:
                logger.error(f"bcftools norm failed: {result_norm.stderr.decode()}")
                return False
        else:
            with open(compressed_vcf, "wb") as f_out:
                result_bgzip = subprocess.run(
                    ["bgzip", "-c"],
                    input=result_view.stdout,
                    stdout=f_out,
                    stderr=subprocess.PIPE
                )
                if result_bgzip.returncode != 0:
                    logger.error(f"bgzip failed: {result_bgzip.stderr.decode()}")
                    return False

        tabix_cmd = ["tabix", "-f", "-p", "vcf", compressed_vcf]
        result_tabix = subprocess.run(
            tabix_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True
        )

        logger.info(f"VCF QC completed: {compressed_vcf}")
        logger.info(f"Tabix index created: {compressed_vcf}.tbi")
        return True

    except subprocess.CalledProcessError as e:
        stderr_output = e.stderr.decode().strip() if e.stderr else "Unknown error"
        logger.error(f"Command failed ({e.cmd}): {stderr_output}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error during VCF QC: {e}")
        return False
