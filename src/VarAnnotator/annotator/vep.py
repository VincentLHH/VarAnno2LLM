from pathlib import Path
import subprocess
from VarAnnotator.utils.logging import setup_logger
import logging

setup_logger()
logger = logging.getLogger(__name__)

def run_vep(
        input_vcf: Path,
        output_vcf: Path,
        vep_cache: Path,
        ref_seq: Path,
        genome: str = "GRCh37",
        enable_plugins: bool = True,
        enable_customs: bool = True,
        enable_pheno: bool = True,
        vep_path: str = None,
        force_overwrite: bool = True,
        output_format: str = "vcf",
        fork: int = 4,
        verbose: bool = True,
) -> bool:
    """ vep annotator"""
    cache_str = ("Uploaded_variation,Location,REF_ALLELE,Allele,ZYG,Consequence,IMPACT,DOMAINS,"
                 "Feature,DISTANCE,EXON,INTRON,SYMBOL,STRAND,HGNC_ID,HGVSc,HGVSp,HGVSg,"
                 "MAX_AF,AF,gnomADg_AF,gnomADg_EAS_AF,Protein_position,Amino_acids,Codons,PUBMED,"
                 "SIFT,PolyPhen,Existing_variation")
    custom_str = ""
    dbnsfp_str = ""
    spliceai_str = ""
    dbscSNV_str = ""
    condel_str = ""
    pheno_str = ""
    input_vcf = Path(input_vcf)
    output_vcf = Path(output_vcf)

    if isinstance(vep_path, str):
        vep_path = vep_path.split()
    if vep_path is None:
        vep_path = ["conda", "run", "-n", "vep", "vep"]

    # validate input and output
    if not input_vcf.exists():
        logger.error(f"Input VCF does not exist: {input_vcf}")
        return False
    if not input_vcf.suffix.lower() in [".vcf", ".gz"]:
        logger.error(f"Input VCF format is not supported: {input_vcf}")
        return False
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    if output_vcf.exists() and not force_overwrite:
        logger.error(f"Output VCF already exists (use force_overwrite): {output_vcf}")
        return False
    # CMD construction
    # plugins = ["spliceAI", "dbNSFP", "dbscSNV", "Condel"] # by default
    # customs = ["clinvar", "intervar"] # by default
    cmd = vep_path
    cmd += [
        "--offline", "--cache", "--merged",
        "--dir_cache", f"{str(vep_cache)}/vep/cache",
        "--dir_plugins", f"{str(vep_cache)}/vep/plugins",
        "--input_file", str(input_vcf),
        "--output_file", str(output_vcf),
        "--keep_csq",
        "--total_length", "--hgvs", "--symbol", "--ccds", "--individual_zyg", "all",
        "--uniprot", "--max_af", "--af", "--af_gnomade", "--af_gnomadg",
        "--sift", "b", "--polyphen", "b",
        "--pick", "--pick_allele", "--pick_order" ,"mane_plus_clinical,mane_select,refseq,canonical,ensembl",
        "--transcript_filter", "stable_id match N[MR]_",
        "--no_escape", "--check_existing", "--exclude_predicted", "--uploaded_allele", "--show_ref_allele",
        "--numbers", "--domains",
        "--format", "vcf",
        f"--{str(output_format)}", "--fa", str(ref_seq),
        "--assembly", str(genome),
        "--fork", str(fork),
    ]
    if force_overwrite:
        cmd.append("--force_overwrite")

    if enable_plugins:
        dbnsfp_str = ("rs_dbSNP,"
                      "REVEL_score,REVEL_rankscore,M-CAP_score,M-CAP_rankscore,M-CAP_pred,GERP++_RS,"
                      "GERP++_RS_rankscore,MVP_score,MVP_rankscore,"
                      "phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,"
                      "phyloP470way_mammalian,phyloP470way_mammalian_rankscore,"
                      "MutationTaster_score,MutationTaster_pred,SIFT_score,SIFT_pred,"
                      "Polyphen2_HDIV_pred,Polyphen2_HVAR_pred")


        spliceai_str = ("SpliceAI_pred_DP_AG,SpliceAI_pred_DP_AL,SpliceAI_pred_DP_DG,SpliceAI_pred_DP_DL,"
                            "SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL")
        if output_format == 'tab':
            spliceai_plugin = (f"SpliceAI,snv={str(vep_cache)}/spliceAI/spliceai_scores.raw.snv.hg19.vcf.gz,indel="
                               f"{str(vep_cache)}/spliceAI/spliceai_scores.raw.indel.hg19.vcf.gz,split_output=1,cutoff=0.5")
        else:
            spliceai_plugin = (f"SpliceAI,snv={str(vep_cache)}/spliceAI/spliceai_scores.raw.snv.hg19.vcf.gz,indel="
                               f"{str(vep_cache)}/spliceAI/spliceai_scores.raw.indel.hg19.vcf.gz,cutoff=0.5")
        dbscSNV_str = "ada_score,rf_score"
        condel_str = "Condel"
        cmd += [
            "--plugin",
            f"dbNSFP,{str(vep_cache)}/dbNSFP/dbNSFP5.2a_grch37.gz,{dbnsfp_str}",
            "--plugin",
            spliceai_plugin,
            "--plugin",
            f"dbscSNV,{str(vep_cache)}/dbscSNV/dbscSNV1.1_GRCh37.txt.gz",
            "--plugin",
            f"Condel,{str(vep_cache)}/vep/plugins/config/Condel/config,b"
        ]
    if enable_customs:
        custom_str = ("ClinVar_CLNSIG,ClinVar_CLNREVSTAT,ClinVar_CLNDN,"
                      "ClinVar_CLNHGVS,InterVar_ACMG,InterVar_Intervar,"
                      "HGMD_HGMD,OMIM")
        cmd += [
            "--custom",
            f"file={str(vep_cache)}/clinvar/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNHGVS",
            "--custom",
            f"file={str(vep_cache)}/intervar/intervar.vcf.gz,short_name=InterVar,format=vcf,type=exact,fields=Intervar%ACMG",
            "--custom",
            f"file={str(vep_cache)}/HGMD/HGMD-2017.vcf.gz,short_name=HGMD,format=vcf,type=exact,fields=HGMD",
            "--custom",
            f"file={str(vep_cache)}/OMIM/omim.hg19.bed.gz,short_name=OMIM,format=bed,type=overlap",
        ]

    if enable_pheno:
        pheno_str = "Pheno_PHENOTYPE,Pheno_ORIG_CPRA"
        cmd += [
            "--custom",
            f"file={str(input_vcf)},short_name=Pheno,format=vcf,type=exact,fields=PHENOTYPE%ORIG_CPRA",
        ]

    field = [cache_str,dbnsfp_str,dbscSNV_str,spliceai_str,condel_str,custom_str,pheno_str]
    field_str = ",".join([f for f in field if f])
    if field_str:
        cmd += [
            "--fields", str(field_str)
        ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False
        )
        if result.returncode == 0:
            logger.info(f"VEP succeeded: {input_vcf} -> {output_vcf}")
            return True
        else:
            logger.error(f"VEP failed with code {result.returncode}")
            logger.error(f"STDERR: {result.stderr}")
            if verbose:
                logger.debug(f"STDOUT: {result.stdout}")
            return False
    except Exception as e:
        logger.error(f"Unexpected error running VEP: {e}")
        return False