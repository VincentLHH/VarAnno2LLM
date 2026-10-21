from pathlib import Path
from types import SimpleNamespace

import pytest

from VarAnnotator.llm.prompt import get_anno, passes_clinical_filters
from VarAnnotator.preprocess.pre import cpra_to_vcf, qc_vcf


def _annotation_record(**overrides):
    record = {
        "Pheno_ORIG_CPRA": "chr1-100-A-C",
        "Pheno_PHENOTYPE": "",
        "Uploaded_variation": "chr1:100:A:C",
        "SYMBOL": "GENE1",
        "Feature": "NM_000001.1",
        "HGVSc": "NM_000001.1:c.10A>C",
        "HGVSp": "NP_000001.1:p.Lys4Thr",
        "Consequence": ["missense_variant"],
        "ZYG": "",
        "ClinVar_CLNSIG": "",
        "HGMD_HGMD": "",
        "InterVar_ACMG": [],
        "InterVar_Intervar": "",
        "SIFT": "",
        "Polyphen2_HDIV_pred": None,
        "Polyphen2_HVAR_pred": None,
        "MutationTaster_pred": "",
        "REVEL_score": "",
        "Condel": "",
        "SpliceAI_pred_DS_AG": "",
        "SpliceAI_pred_DS_AL": "",
        "SpliceAI_pred_DS_DG": "",
        "SpliceAI_pred_DS_DL": "",
        "rf_score": "",
        "ada_score": "",
        "phyloP100way_vertebrate": "",
        "phyloP470way_mammalian": "",
        "GERP++_RS": "",
        "OMIM": "",
        "AF": "",
        "gnomADg_AF": "",
        "gnomADg_EAS_AF": "",
        "MAX_AF": "",
    }
    record.update(overrides)
    return record


def test_spliceai_max_uses_non_ag_delta_when_ag_is_missing():
    """Using AG presence as a gate for every SpliceAI branch loses valid scores."""
    annotation = _annotation_record(
        SpliceAI_pred_DS_AG="",
        SpliceAI_pred_DS_AL="0.12",
        SpliceAI_pred_DS_DG="0.44",
        SpliceAI_pred_DS_DL="0.01",
    )

    assert get_anno(annotation)["SpliceAI"] == pytest.approx(0.44)


def test_cpra_without_zygosity_defaults_to_heterozygous(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """The documented empty-zygosity default must be reflected in generated GT."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "input.tsv").write_text("chr1-100-A-C\tphenotype\t\n")
    monkeypatch.setattr(
        "VarAnnotator.preprocess.pre.subprocess.run",
        lambda *args, **kwargs: SimpleNamespace(returncode=0, stdout=b"", stderr=b""),
    )

    assert cpra_to_vcf(tmp_path / "input.tsv", Path("output.vcf"))
    generated = Path("unsorted_output.vcf").read_text()
    assert "\tGT:AD:DP\t0/1:50,50:100" in generated


def test_qc_vcf_splits_multiallelic_records(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """Without -m -any, multiallelic records are not decomposed."""
    input_vcf = tmp_path / "input.vcf"
    reference = tmp_path / "reference.fa"
    output_vcf = tmp_path / "output"
    input_vcf.write_text("##fileformat=VCFv4.2\n")
    reference.write_text(">chr1\nA\n")
    calls = []

    def fake_run(command, **kwargs):
        calls.append(command)
        return SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr("VarAnnotator.preprocess.pre.subprocess.run", fake_run)

    assert qc_vcf(input_vcf, output_vcf, reference)
    norm_command = next(command for command in calls if command[:2] == ["bcftools", "norm"])
    assert norm_command[2:6] == ["-c", "s", "-m", "-any"]


@pytest.mark.parametrize(
    ("consequence", "expected"),
    [
        (["missense_variant"], True),
        (["stop_gained"], True),
        (["splice_donor_variant"], True),
        (["synonymous_variant"], False),
        (["intron_variant", "splice_acceptor_variant"], True),
    ],
)
def test_consequence_filter_keeps_protein_or_splice_altering_terms(
    consequence, expected
):
    """The clinical filter must not treat all VEP consequences as equivalent."""
    result = passes_clinical_filters(
        cpra="chr1-100-A-C",
        consequence=consequence,
        associated_diseases="",
        consequence_filter=True,
        omim_filter=False,
        max_variant_length=None,
    )
    assert result is expected


def test_clinical_filter_rejects_very_long_or_non_omim_records():
    """The <1 kb and OMIM filters described in the pipeline must be enforceable."""
    long_cpra = "chr1-100-" + ("A" * 1000) + "-C"
    assert not passes_clinical_filters(
        cpra=long_cpra,
        consequence=["missense_variant"],
        associated_diseases='{"description":"disease","genemap":[]}',
        consequence_filter=True,
        omim_filter=True,
        max_variant_length=1000,
    )
    assert not passes_clinical_filters(
        cpra="chr1-100-A-C",
        consequence=["missense_variant"],
        associated_diseases="",
        consequence_filter=True,
        omim_filter=True,
        max_variant_length=1000,
    )
