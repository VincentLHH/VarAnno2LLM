import argparse
from VarAnnotator.utils.logging import setup_logger
import logging
from VarAnnotator.annotator.vep import run_vep
from VarAnnotator.preprocess.pre import cpra_to_vcf, qc_vcf
from VarAnnotator.llm.prompt import get_jsonl
from humanfriendly.terminal import output

setup_logger()
logger = logging.getLogger(__name__)

def run_annotator(
        args: argparse.ArgumentParser
) -> bool:


    if not args.input.exists():
        logger.error(f"Input file not exist: {args.input}")
        return False

    if not args.vep_cache.exists():
        logger.error(f"VEP cache path not exist: {args.vep_cache}")
        return False

    if not args.ref_seq.exists():
        logger.error(f"Ref file not exist: {args.ref_seq}")
        return False

    output_path = args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Start annotating: {args.input} → {output_path}")

    success = run_vep(
        input_vcf=args.input,
        output_vcf=output_path,
        vep_cache=args.vep_cache,
        ref_seq=args.ref_seq,
        genome=args.genome,
        vep_path=args.vep_path,
        force_overwrite=not args.no_overwrite,
        enable_plugins=not args.no_plugins,
        enable_customs=not args.no_customs,
        fork=args.fork,
        verbose=args.verbose,
        output_format=args.output_format,
    )

    if success:
        logger.info(f"Done for annotation: {output_path}")
        return True
    else:
        logger.error(f"Failed to annotate, check the log plz")
        return False

def run_preprocessor(
        args: argparse.ArgumentParser,
) -> bool:
    logger.info(f"Start qc and normalization: {args.input} → {args.output}")
    if not args.input.exists():
        logger.error(f"Input file not exist: {args.input}")
        return False
    if args.cpra:
        success_cpra = cpra_to_vcf(
            cpra_input=args.input,
            output_vcf=args.output,
        )
        if not success_cpra:
            logger.error("Failed to convert CPRA to VCF")
            return False
    raw_vcf = args.output.with_suffix(".vcf")
    compressed_vcf = str(raw_vcf) + ".gz"
    if not args.no_qc:
        if not args.ref_seq.exists():
            logger.error(f"Ref file not exist: {args.ref_seq}")
            return False
        success_qc = qc_vcf(
            input_vcf=args.input if not args.cpra else compressed_vcf,
            output_vcf=args.output,
            ref_dir=args.ref_seq,
            min_qual=args.min_qual,
        )
        if not success_qc:
            logger.error("Failed to QC VCF")
            return False
    logger.info(f"Done for preprocessing: {args.output}" + ".gz")
    return True

def run_prompter(
        args: argparse.ArgumentParser,
) -> bool:
    logger.info(f"Start prompting: {args.input} → {args.output}")
    if not args.input.exists():
        logger.error(f"Input file not exist: {args.input}")
        return False
    if not args.gene_list.exists():
        logger.error(f"Gene list file not found: {args.gene_list}")
        return False
    success = get_jsonl(
        input = args.input,
        output = args.output,
        gene_set_file = args.gene_list,
        white_list_file = args.white_list if args.white_list else None,
        pheno_file = args.pheno_file if args.pheno_file else None,
        label = args.label,
        af_filter = args.af_filter,
        use_white_list = not args.no_whitelist,
        consequence_filter = args.consequence_filter,
        omim_filter = args.omim_filter,
        max_variant_length = args.max_variant_length,
    )
    if success:
        logger.info(f"Done for prompt: {args.output}/{args.input.name}.jsonl")
        return True
    else:
        logger.error(f"Failed to prompt, check the log plz")
        return False
