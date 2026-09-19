"""Shared fixtures.

The suite runs evomining-py over the EvoMining "los17" example genomes and compares the
result against committed snapshots (see ``golden.py``).  There is no second
implementation to check against -- evomining-py's own current output is the baseline, so
these tests answer "did anything change?", not "is this biologically right?".

Two tiers, by what they need:

* **fast** -- ``analyze`` replayed from committed BLAST output and trimmed metadata
  (~130 KB). No genomes, no external tools, runs in seconds.
* **`genomes` / `tools`** -- the stages that read the 233 MB genome set or shell out to
  BLAST, MUSCLE and FastTree. These skip with instructions when their inputs are absent.
"""

from __future__ import annotations

import shutil
from types import SimpleNamespace

import pytest

from evomining import analysis, genomedb, start, trees
from evomining.external import Tools

DATA = __import__("pathlib").Path(__file__).resolve().parent / "data"

GENOMES_DIR = DATA / "los17_genomes"
ENZYME_DB = DATA / "enzymes_db.faa"
BLAST_DIR = DATA / "blast"
METADATA_DIR = DATA / "metadata"

#: MIBiG protein set, downloaded separately (~31 MB, not committed). Snapshots were taken
#: with version 4.0; a different release will legitimately move the MIBiG-derived results.
MIBIG_GLOB = "mibig_prot_seqs_*.fasta"

#: How many genomes the example set holds. `accessions.txt` lists them; one of the
#: original 18 los17 was not included in the ORIGINAL evomining test analysis
#: 3 others were renamed since. The correct accessions are in the accessions.txt and
#: are in total 17
N_GENOMES = 17


def pytest_addoption(parser):
    parser.addoption(
        "--update-golden", action="store_true", default=False,
        help="rewrite the snapshots under tests/golden/ instead of comparing to them")


def pytest_configure(config):
    config.addinivalue_line("markers", "genomes: needs tests/data/los17_genomes/")
    config.addinivalue_line("markers", "mibig: needs the MIBiG protein FASTA")


def require_tools(*names):
    """Skip unless every named external tool resolves on $PATH."""
    tools = Tools()
    missing = [n for n in names if not tools.available(n)]
    if missing:
        pytest.skip(f"external tool(s) not on $PATH: {', '.join(missing)}")


# ---------------------------------------------------------------------------
#  Inputs
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def genomes_dir():
    """The GenBank example set, or skip explaining how to fetch it."""
    if not GENOMES_DIR.is_dir() or not any(GENOMES_DIR.glob("*.gbff")):
        pytest.skip(
            f"example genomes not found at {GENOMES_DIR}.\n"
            f"Fetch them with the accessions in {DATA / 'accessions.txt'}:\n"
            "  datasets download genome accession --inputfile tests/data/accessions.txt "
            "--include gbff --filename los17_gbff.zip\n"
            "then place one <accession>.gbff per genome in tests/data/los17_genomes/.")
    return GENOMES_DIR


@pytest.fixture(scope="session")
def mibig_fasta():
    """The MIBiG protein FASTA, or skip explaining how to get it.

    Downloaded separately rather than committed: it is ~31 MB, it is a versioned public
    reference set, and pinning a trimmed copy would quietly change BLAST e-values (which
    scale with database size) relative to a real run.
    """
    found = sorted(DATA.glob(MIBIG_GLOB))
    if not found:
        pytest.skip(
            f"MIBiG protein FASTA not found at {DATA}/{MIBIG_GLOB}.\n"
            "Download it from https://mibig.secondarymetabolites.org/download\n"
            "(the 'MIBiG proteins' FASTA) and place it in tests/data/.\n"
            "Snapshots were generated with MIBiG 4.0.")
    return found[-1]


# ---------------------------------------------------------------------------
#  The fast tier: analyze, replayed from committed BLAST
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def analysis_out(tmp_path_factory):
    """Run `analyze` from the committed BLAST output and trimmed metadata.

    The metadata tables are trimmed to just the proteins the BLAST files mention. That is
    verified to produce byte-identical output to running against the full 86k-protein
    tables, because `analyze` only ever looks up proteins it has a hit for -- which is
    what lets this tier run without the genome set.

    `--mibig` is deliberately omitted: it is the only part of `analyze` that shells out
    to BLAST, and leaving it off keeps this tier tool-free.
    """
    work = tmp_path_factory.mktemp("analyze")
    blast = work / "blast"
    blast.mkdir()
    for name in ("thereand.blast", "backagain.blast"):
        shutil.copy(BLAST_DIR / name, blast / name)

    out = work / "evomining_analysis"
    analysis.run(SimpleNamespace(
        blast_dir=str(blast),
        central_db=str(ENZYME_DB),
        genome_names=str(METADATA_DIR / "genome_names.tsv"),
        genome_functions=str(METADATA_DIR / "genome_functions.tsv"),
        genomes_fasta=str(ENZYME_DB),   # unused without --mibig; must merely exist
        antismash=None, mibig=None,
        output_dir=str(out), threads=2, sd=1.0, ccm="multi",
    ), tools=Tools())
    return out


# ---------------------------------------------------------------------------
#  The slow tier: real runs over the genome set
# ---------------------------------------------------------------------------

def _require_pipeline_tools():
    require_tools("makeblastdb", "blastp", "muscle", "trimal",
                  "nw_clade", "nw_labels", "nw_reroot", "nw_rename", "nw_display")
    tools = Tools()
    if not (tools.available("FastTree") or tools.available("fasttree")):
        pytest.skip("FastTree not on $PATH")
    return tools


@pytest.fixture(scope="session")
def pipeline(genomes_dir, tmp_path_factory):
    """genome-db -> start -> analyze, over the real genome set. No MIBiG.

    Session-scoped and shared with the MIBiG tests, so the ~40 s of genome parsing and
    BLAST happens once per run rather than once per module.
    """
    tools = _require_pipeline_tools()
    work = tmp_path_factory.mktemp("pipeline")
    db = work / "evomining_db"

    genomedb.run(SimpleNamespace(
        input_dir=str(genomes_dir), output_dir=str(db),
        lists=None, list_dir=None, names=None, keep_pseudogenes=False))

    start.run(SimpleNamespace(
        genomes=str(db / "GENOMES.fasta"), central_db=str(ENZYME_DB),
        output_dir=str(work / "evomining_out"),
        threads=4, evalue=1e-4, chunk_size=100_000), tools=tools)
    blast_dir = next((work / "evomining_out").glob("*/blast"))

    analysis_dir = work / "evomining_analysis"
    analysis.run(SimpleNamespace(
        blast_dir=str(blast_dir), central_db=str(ENZYME_DB),
        genome_names=str(db / "genome_names.tsv"),
        genome_functions=str(db / "genome_functions.tsv"),
        genomes_fasta=str(db / "GENOMES.fasta"),
        antismash=None, mibig=None, output_dir=str(analysis_dir),
        threads=4, sd=1.0, ccm="multi"), tools=tools)

    return SimpleNamespace(work=work, db=db, blast=blast_dir,
                           analysis=analysis_dir, tools=tools)


@pytest.fixture(scope="session")
def mibig_run(pipeline, mibig_fasta):
    """analyze --mibig then trees --mibig, on top of `pipeline`.

    This is the only path that exercises `find_evomining_predictions`: without MIBiG
    reference leaves the algorithm has no seeds to walk out from and returns nothing, so
    the green "EvoMining prediction" class is unreachable.
    """
    work = pipeline.work
    analysis_dir = work / "analysis_mibig"
    analysis.run(SimpleNamespace(
        blast_dir=str(pipeline.blast), central_db=str(ENZYME_DB),
        genome_names=str(pipeline.db / "genome_names.tsv"),
        genome_functions=str(pipeline.db / "genome_functions.tsv"),
        genomes_fasta=str(pipeline.db / "GENOMES.fasta"),
        antismash=None, mibig=str(mibig_fasta), output_dir=str(analysis_dir),
        threads=4, sd=1.0, ccm="multi"), tools=pipeline.tools)

    trees_dir = work / "trees_mibig"
    trees.run(SimpleNamespace(
        classified=str(analysis_dir / "classified_proteins.tsv"),
        central_db=str(ENZYME_DB), genomes_fasta=str(pipeline.db / "GENOMES.fasta"),
        mibig=str(mibig_fasta),
        mibig_blast=str(analysis_dir / "expanded_vs_mibig.blast"),
        families=None, itol=True, microreact=True, skip_existing=False,
        max_seqs=None, threads=4, output_dir=str(trees_dir)), tools=pipeline.tools)

    return SimpleNamespace(analysis=analysis_dir, trees=trees_dir, mibig=mibig_fasta)


@pytest.fixture(scope="session")
def genome_names():
    return analysis.load_genome_names(METADATA_DIR / "genome_names.tsv")


@pytest.fixture(scope="session")
def forward_hits():
    return analysis.parse_blast_results(str(BLAST_DIR / "thereand.blast"))


@pytest.fixture(scope="session")
def reverse_hits():
    return analysis.parse_blast_results(str(BLAST_DIR / "backagain.blast"))
