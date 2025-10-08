from VarAnnotator.utils.logging import setup_logger
from VarAnnotator.utils.io import get_anno_dict, get_gene_set, get_pheno_info, get_white_list
import logging
from pathlib import Path
from tqdm import tqdm
import pandas as pd
import json

setup_logger()
logger = logging.getLogger(__name__)

def get_anno(
        anno: dict,
) -> dict:
    snv_list = {}

    # Base Info
    tmp = anno["Pheno_ORIG_CPRA"] if anno["Pheno_ORIG_CPRA"] else ""
    cpra1 = anno["Pheno_ORIG_CPRA"][0] if isinstance(anno["Pheno_ORIG_CPRA"], (list, tuple)) else tmp
    cpr = anno["Uploaded_variation"].split("/")[0].replace(":", "-").replace("_", "-")
    a = anno["Uploaded_variation"].split("/")[-1].replace(":", "-").replace("_", "-")
    cpra2 = "-".join([cpr, a])
    snv_list["cpra"] = cpra1 if cpra1 else cpra2
    snv_list["phenotype"] = anno["Pheno_PHENOTYPE"].replace('_', ',') if anno["Pheno_PHENOTYPE"] else ""
    snv_list["gene symbol"] = anno["SYMBOL"] if anno["SYMBOL"] else ""
    snv_list["transcript"] = anno["Feature"] if anno["Feature"] else ""
    snv_list["dna change"] = anno["HGVSc"].split(":")[1] if anno["HGVSc"] and ":" in anno["HGVSc"] else ""
    snv_list["protein change"] = anno["HGVSp"].split(":")[1] if anno["HGVSp"] and ":" in anno["HGVSp"] else ""
    snv_list["functional consequences"] = anno["Consequence"] if anno["Consequence"] else ""
    snv_list["zygosity"] = anno["ZYG"].split(":")[1].capitalize() if anno["ZYG"] else ""

    # Database
    snv_list["Clinvar database"] = anno["ClinVar_CLNSIG"] if anno["ClinVar_CLNSIG"] else ""
    snv_list["HGMD database"] = anno["HGMD_HGMD"] if anno["HGMD_HGMD"] else ""
    snv_list["ACMG evidence items"] = anno["InterVar_ACMG"] if anno["InterVar_ACMG"] else ""
    snv_list["pathogenicity assessment"] = anno["InterVar_Intervar"] if anno["InterVar_Intervar"] else ""

    # predictor
    snv_list["SIFT"] = anno["SIFT"].split("(")[0] if anno["SIFT"] and "(" in anno["SIFT"] else ""
    Polyphen_amp = {
        "D": "probably damaging",
        "P": "possibly damaging",
        "B": "benign",
    }
    snv_list["Polyphen2-HDIV"] = Polyphen_amp.get(anno["Polyphen2_HDIV_pred"][0],"") if anno["Polyphen2_HDIV_pred"] else ""
    snv_list["Polyphen2-HVAR"] = Polyphen_amp.get(anno["Polyphen2_HVAR_pred"][0],"") if anno["Polyphen2_HVAR_pred"] else ""
    MutationTaster_map = {
        "A": "disease_causing_automatic",
        "D": "disease_causing",
        "N": "polymorphism",
        "P": "polymorphism_automatic"
    }
    MT2 = anno.get("MutationTaster_pred", "")
    snv_list["MutationTaster2"] = MutationTaster_map.get(MT2[0], "") if isinstance(MT2, str) and len(MT2) > 0 else ""
    snv_list["Revel"] = anno["REVEL_score"][0] if anno["REVEL_score"] else ""
    snv_list["Condel"] = anno["Condel"].split('(')[1].split(')')[0] if anno["Condel"] else ""
    splice_scores = [
        float(anno["SpliceAI_pred_DS_AG"] if anno["SpliceAI_pred_DS_AG"] else 0),
        float(anno["SpliceAI_pred_DS_AL"] if anno["SpliceAI_pred_DS_AG"] else 0),
        float(anno["SpliceAI_pred_DS_DG"] if anno["SpliceAI_pred_DS_AG"] else 0),
        float(anno["SpliceAI_pred_DS_DL"] if anno["SpliceAI_pred_DS_AG"] else 0),
    ]
    splice_scores = [s for s in splice_scores if not pd.isna(s) and s >= 0]
    snv_list["SpliceAI"] = max(splice_scores) if splice_scores else ""
    snv_list["dbscSnv_RF"] = anno["rf_score"] if anno["rf_score"] else ""
    snv_list["dbscSnv_ADA"] = anno["ada_score"] if anno["ada_score"] else ""
    snv_list["PhyloP Vertebrates"] = anno["phyloP100way_vertebrate"] if anno["phyloP100way_vertebrate"] else ""
    snv_list["PhyloP Mammals"] = anno["phyloP470way_mammalian"] if anno["phyloP470way_mammalian"] else ""
    snv_list["GERP++"] = anno["GERP++_RS"] if anno["GERP++_RS"] else ""

    # disease
    if anno["OMIM"]:
        des = anno["OMIM"].split("|")[2] if anno["OMIM"].split("|")[2] and anno["OMIM"].split("|")[2] != "-" else ""
        gmp = anno["OMIM"].split("|")[3] if anno["OMIM"].split("|")[3] and anno["OMIM"].split("|")[3] != "-" else ""
        genemap = ""
        if gmp:
            try:
                genemap = json.loads(gmp.replace("'", '"'))
            except Exception as e:
                print(e)

        snv_list["associated diseases"] = str({"description": des, "genemap": genemap})
        snv_list["inheritance"]= ",".join(set([item.get('Inheritance', '-') for item in genemap]))
    else:
        snv_list["associated diseases"] = ""
        snv_list["inheritance"] = ""

    # AF
    snv_list["1000G"] = anno["AF"] if anno["AF"] else ""
    snv_list["gnomAD"] = anno["gnomADg_AF"] if anno["gnomADg_AF"] else ""
    snv_list["gnomAD-EAS"] = anno["gnomADg_EAS_AF"] if anno["gnomADg_EAS_AF"] else ""
    snv_list["max_af"] = anno["MAX_AF"] if anno["MAX_AF"] else ""

    # ACMG
    snv_list["ACMG evidence items"] = ','.join(anno["InterVar_ACMG"]) if anno["InterVar_ACMG"] else ""
    snv_list["pathogenicity assessment"] = anno["InterVar_Intervar"] if anno["InterVar_Intervar"] else ""

    return snv_list

def get_prompt(
        filtered_dict: dict,
        phenotype: str,
        gene_set: set,
        label: str = "test",
) -> dict:
    labels = ["test", "positive", "negative"]
    if label not in labels:
        logger.error(f"{label} is not a valid label in {labels}")
        return {}

    # IEI
    # phenotype = "Primary immunodeficiency disease, inborn errors of immunity, " + phenotype
    phenotype = "Inherited retinal disease, " + phenotype

    # IEI
    # system_context = ("Please play the role of a clinical genetic analyst that involves analyzing the Next-Generation "
    #                   "Sequencing (NGS) data of a patient with suspected inborn errors of immunity (IEI) or primary immunodeficiency disease (PID), "
    #                   "considering their clinical phenotype and annotation information to identify genetic variations that may "
    #                   "contribute to the observed clinical phenotype.")

    # IRD
    system_context = ("Please play the role of a clinical genetic analyst that involves analyzing the Next-Generation "
                      "Sequencing (NGS) data of a patient with suspected inherited retinal disease (IRD), "
                      "considering their clinical phenotype and annotation information to identify genetic variations that may "
                      "contribute to the observed clinical phenotype.")

    quest_begin = ("What are your responsibilities as a medical expert in the diagnostic process and what are the keys "
                   "to succeeding in making accurate diagnoses?")
    response_begin = (f"As a seasoned diagnostic expert, I focus on biological pathogenicity, inheritance mode, "
                      f"and clinical phenotype to identify which variations truly cause disease or phenotypic manifestations."
                      f"\nThe basic information of this variant is as follows: "
                      f"\nGene symbol: {filtered_dict['gene symbol']} \nTranscript: {filtered_dict['transcript']} "
                      f"\nDNA change: {filtered_dict['dna change']} \nProtein change: {filtered_dict['protein change']} ")
    quest_thought = "What is your plan at current step?"
    response1_1 = (f"Element 1, I will consider the biological hazard of the mutation according to the database "
                   f"query results, predictions from software and functional changes.")
    quest_meaning = "How do you interpret the information at current step?"
    response1_2 = (
        "Firstly, the biological significance from ClinVar ranges from strong to weak: Pathogenic, Likely pathogenic, "
        "Likely benign, Benign/Likely_benign, Benign. Special cases: Uncertain_significance, "
        "Conflicting_interpretations_of_pathogenicity. In HGMD, \"DM\" or \"DM?\" indicates strong "
        "significance; \"DFP\", \"DP\", or \"FP\" indicates weak significance.\nSecondly, for software predictions ("
        "SIFT, Polyphen2-HDIV, etc.), \"deleterious\" results indicate strong significance; \"polymorphism\" "
        "indicates weak. Conservation tools (PhyloP, GERP++) showing \"conserved\" mark strong significance; "
        "\"non-conserved\" marks weak.\nFunctional changes also impact significance, including categories such "
        "as synonymous (no change in protein), missense (amino acid change), nonsense (premature stop codon), "
        "frameshift (disrupting protein structure) and so on.\nThirdly, other specific values and meanings include: "
        "\n -REVEL, Range 0 to 1 (Higher scores indicate stronger pathogenicity). "
        "\n -Condel, Range 0 to 1 (Higher values mean higher deleterious potential). "
        "\n -SpliceAI, Range 0 to 1 (Higher scores indicate stronger splicing impact).\n -dbscSNV_RF/ADA: Range 0 "
        "to 1 (Higher scores indicate a stronger splicing impact).\n -PhyloP Vertebrates/Mammals, "
        "Range negative to positive (Higher values mean higher conservation; Negative values indicate acceleration). "
        "\n -GERP++, Range negative to positive (Higher values indicate stronger evolutionary "
        "constraint).")
    quest_value = "What are the values of these terms?"
    response1_3 = (f"For element 1, the information of the variant is as follows: "
                   f"\nClinvar database: {filtered_dict['Clinvar database']} "
                   f"\nHGMD database: {filtered_dict['HGMD database']} "
                   f"\nACMG evidence items: {filtered_dict['ACMG evidence items']} "
                   f"\nDatabase assessment: {filtered_dict['pathogenicity assessment']} "
                   f"\nFunctional consequences: {filtered_dict['functional consequences']} "
                   f"\nSIFT: {filtered_dict['SIFT']} "
                   f"\nPolyphen2-HDIV: {filtered_dict['Polyphen2-HDIV']} "
                   f"\nPolyphen2-HVAR: {filtered_dict['Polyphen2-HVAR']} "
                   f"\nMutationTaster2: {filtered_dict['MutationTaster2']} "
                   f"\nRevel: {filtered_dict['Revel']} "
                   f"\nCondel: {filtered_dict['Condel']} "
                   f"\nSpliceAI: {filtered_dict['SpliceAI']} "
                   f"\ndbscSnv_RF: {filtered_dict['dbscSnv_RF']} "
                   f"\ndbscSnv_ADA: {filtered_dict['dbscSnv_ADA']} "
                   f"\nPhyloP Vertebrates: {filtered_dict['PhyloP Vertebrates']} "
                   f"\nPhyloP Mammals: {filtered_dict['PhyloP Mammals']} "
                   f"\nGERP++: {filtered_dict['GERP++']} ")
    response2_1 = (f"Element 2, I will consider the inheritance mode of the mutation according to the "
                   f"zygosity and allele frequencies from databases of normal individuals.")
    response2_2 = ("Firstly, I classify the zygosity of a mutation into homozygous (Hom), heterozygous (Het), "
                   "and hemizygous (Hemi). The inheritance mode of a mutation can be classified into autosomal dominant ("
                   "AD), autosomal recessive (AR), X-linked dominant (XLD), X-linked recessive (XLR), and X-linked ("
                   "XL).\nSecondly, I assess whether the patient's zygosity conforms to genetic inheritance patterns: If "
                   "zygosity is \"Hom\" and the genetic mode is \"AR\", I consider it consistent; if zygosity is \"Het\" "
                   "and the genetic mode is \"AD\", I consider it consistent (note rare cases of \"Hom\" under \"AD\" "
                   "indicating more severe phenotypes); if zygosity is \"Hemi\" or \"Het\" and the genetic mode is "
                   "\"XL\", I consider it consistent.\nThirdly, I refer to database query results for variant carriage in "
                   "normal individuals: In databases such as 1000G, gnomAD, and gnomAD-EAS, "
                   "if the result is not empty, normal individuals carry the variant. If the allele frequency (AF) is "
                   "less than 0.05, I consider the variant rare; if AF is less than 0.01 or the result is empty, "
                   "I consider the variant extremely rare. If normal individuals carry the variant, I also consider "
                   "whether it leads to subtle phenotypes, late onset, or incomplete penetrance, or if the patient's "
                   "zygosity does not match the expected genetic mode for the gene.")
    response2_3 = (f"For element 2, the information of the variant is as follows: "
                   f"\nZygosity: {filtered_dict['zygosity']} "
                   f"\nInheritance mode: {filtered_dict['inheritance']} "
                   f"\n1000G: {filtered_dict['1000G']} "
                   f"\ngnomAD: {filtered_dict['gnomAD']} "
                   f"\ngnomAD-EAS: {filtered_dict['gnomAD-EAS']} ")

    # IEI
    # response3_1 = (f"Element 3, I will consider whether the gene corresponding to the variant occurs in "
    #                f"the IEI/PID gene list provided by IUIS to support the diagnosis.")
    # response3_2 = ("I use the list of pathogenic genes for IEI/PID provided by IUIS based on prior research. If a "
    #                "gene is on this list, I consider the variant more likely to be pathogenic, causing the "
    #                "clinical phenotype.")
    # answer_rag = "yes" if filtered_dict['gene symbol'] in gene_set else "no"
    # response3_3 = (f"For element 3, the information of variant {filtered_dict['cpra']} is as follows: "
    #                f"\nIf \"{filtered_dict['gene symbol']}\" is on the IEI/PID gene list: {answer_rag}")

    # IRD
    response3_1 = (f"Element 3, I will consider whether the gene corresponding to the variant occurs in "
                   f"the IRD gene list provided by RetNet to support the diagnosis.")
    response3_2 = ("I use the list of pathogenic genes for IRD provided by RetNet based on prior research. If a "
                   "gene is on this list, I consider the variant more likely to be pathogenic, causing the "
                   "clinical phenotype.")
    answer_rag = "on" if filtered_dict['gene symbol'] in gene_set else "not on"
    response3_3 = (f"For element 3, the information of the variant is as follows: "
                   f"\nGene \"{filtered_dict['gene symbol']}\" is {answer_rag} the IRD gene list.")

    response4_1 = (f"Element 4, I will consider the clinical phenotype and determine "
                   "whether the information of each disease in the disease list from OMIM is related to clinical phenotypes "
                   "and diagnoses, including disease name, disease details, disease inheritance, and disease phenotype.")
    # IEI
    # response4_2 = ("Firstly, I extract the disease information related to the variant from OMIM, including the "
    #                "disease name, MIM number, disease description, inheritance mode, and other relevant "
    #                "details.\nSecondly, I compare the extracted disease information directly with the patient's "
    #                "clinical phenotype, checking for matching phenotypes.\nThirdly, "
    #                "I further investigate the potential biological associations between the disease information and the "
    #                "patient's phenotype to avoid missing connections due to inconsistent naming but clinically "
    #                "related phenotypes.\nAdditionally, considering the heterogeneity and subtlety of IEI "
    #                "phenotypes, I also consider the potential association between the variant and phenotype in "
    #                "cases where there is no direct match.")

    # IRD
    response4_2 = ("Firstly, I extract the disease information related to the variant from OMIM, including the "
                   "disease name, MIM number, disease description, inheritance mode, and other relevant "
                   "details.\nSecondly, I compare the extracted disease information directly with the patient's "
                   "clinical phenotype, checking for matching phenotypes.\nThirdly, "
                   "I further investigate the potential biological associations between the disease information and the "
                   "patient's phenotype to avoid missing connections due to inconsistent naming but clinically "
                   "related phenotypes.\nAdditionally, considering the heterogeneity and subtlety of IRD "
                   "phenotypes, I also consider the potential association between the variant and phenotype in "
                   "cases where there is no direct match.")
    # IEI
    # response4_3 = (f"For element 4, the information of variant {filtered_dict['cpra']} is as follows: "
    #                f"\nDisease JSON: {filtered_dict['associated diseases']} "
    #                f"\nClinical phenotype: Primary immunodeficiency disease, inborn errors of immunity, {phenotype}")

    # IRD
    response4_3 = (f"For element 4, the information of the variant is as follows: "
                   f"\nDisease JSON: {filtered_dict['associated diseases']} "
                   f"\nClinical phenotype: {phenotype}")

    response5_1 = "Lastly, I will integrate all the data and previous analyses to make the final diagnosis."
    response5_2 = (
        f"I will only reply a single word \"positive\" if the variant causes any of the phenotypes "
        f"\"{phenotype}\", and \"negative\" if it does not. I will not say any other words.")
    quest_final = "What is your final diagnosis?"
    anno = json.dumps(filtered_dict)
    tag = str(anno) if label == "test" else label
    pos_dict = {
        "messages": [
            {
                "role": "system",
                "content": str(system_context),
            },
            {
                "role": "user",
                "content": str(quest_begin),
            },
            {
                "role": "assistant",
                "content": str(response_begin),
            },
            {
                "role": "user",
                "content": str(quest_thought),
            },
            {
                "role": "assistant",
                "content": str(response1_1),
            },
            {
                "role": "user",
                "content": str(quest_meaning),
            },
            {
                "role": "assistant",
                "content": str(response1_2),
            },
            {
                "role": "user",
                "content": str(quest_value),
            },
            {
                "role": "assistant",
                "content": str(response1_3),
            },
            {
                "role": "user",
                "content": str(quest_thought),
            },
            {
                "role": "assistant",
                "content": str(response2_1),
            },
            {
                "role": "user",
                "content": str(quest_meaning),
            },
            {
                "role": "assistant",
                "content": str(response2_2),
            },
            {
                "role": "user",
                "content": str(quest_value),
            },
            {
                "role": "assistant",
                "content": str(response2_3),
            },
            {
                "role": "user",
                "content": str(quest_thought),
            },
            {
                "role": "assistant",
                "content": str(response3_1),
            },
            {
                "role": "user",
                "content": str(quest_meaning),
            },
            {
                "role": "assistant",
                "content": str(response3_2),
            },
            {
                "role": "user",
                "content": str(quest_value),
            },
            {
                "role": "assistant",
                "content": str(response3_3),
            },
            {
                "role": "user",
                "content": str(quest_thought),
            },
            {
                "role": "assistant",
                "content": str(response4_1),
            },
            {
                "role": "user",
                "content": str(quest_meaning),
            },
            {
                "role": "assistant",
                "content": str(response4_2),
            },
            {
                "role": "user",
                "content": str(quest_value),
            },
            {
                "role": "assistant",
                "content": str(response4_3),
            },
            {
                "role": "user",
                "content": str(quest_thought),
            },
            {
                "role": "assistant",
                "content": str(response5_1),
            },
            {
                "role": "user",
                "content": str(quest_meaning),
            },
            {
                "role": "assistant",
                "content": str(response5_2),
            },
            {
                "role": "user",
                "content": str(quest_final),
            },
            {
                "role": "assistant",
                "content": tag,
            },
        ],
    }
    #print(pos_dict)
    return pos_dict

def get_zyg(
        inh: list,
        zyg_ori: str,
) -> str:
    zyg = zyg_ori
    if not inh:
        return zyg
    else:
        if "AD" in inh:
            zyg = "Het"
        elif "DN" in inh:
            zyg = "Het"
        elif "AR" in inh:
            zyg = "Hom"
        elif "XLD" in inh:
            zyg = "Het"
        elif "XLR" in inh:
            zyg = "Hemi"
        elif "XL" in inh:
            zyg = "Hemi"
        else:
            zyg = zyg_ori
    return zyg

def get_jsonl(
        input: Path,
        output: Path,
        gene_set_file: Path,
        white_list_file: Path = None,
        pheno_file: Path = None,
        label: str = "test",
        af_filter: bool = False,
        use_white_list: bool = True,
) -> bool:
    gene_set = get_gene_set(gene_set_file)
    white_list = get_white_list(white_list_file) if white_list_file else None
    pheno_file = get_pheno_info(pheno_file) if pheno_file else None
    sample_name = input.name
    output_file = output.joinpath(sample_name + ".jsonl")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(output_file, 'w', encoding='utf-8') as f_out:
            for record in tqdm(get_anno_dict(input), desc="Prompting"):
                anno = get_anno(record)
                if af_filter and use_white_list and white_list:
                    max_af = float(anno["max_af"]) if anno["max_af"] else 0
                    if max_af > 0.01 and anno["cpra"] not in white_list:
                        continue
                elif af_filter and not use_white_list:
                    max_af = float(anno["max_af"]) if anno["max_af"] else 0
                    if max_af > 0.01:
                        continue

                if label != "test":
                    inh_list = [x.strip() for x in anno["inheritance"].split(",")] if anno["inheritance"] else []
                    anno["zygosity"] = get_zyg(inh_list, anno["zygosity"])
                # cpra = anno["cpra"]
                pheno = anno["phenotype"]
                if pheno_file:
                    phenotype = pheno_file[f"{sample_name}"]
                else:
                    phenotype = pheno
                prompt = get_prompt(anno, phenotype, gene_set, label)
                f_out.write(json.dumps(prompt, ensure_ascii=False) + '\n')
        return True
    except Exception as e:
        logger.error(e)
        return False