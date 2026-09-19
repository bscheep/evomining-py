"""genome-db -> start -> analyze, over the real genome set.

Needs the genomes and the toolchain.  Its most important job is the middle: it re-derives
the BLAST output and metadata that ``test_analyze.py`` replays from committed fixtures.
Without this, those fixtures could drift out of step with the code and the fast tier would
keep passing against a stale baseline.
"""

from __future__ import annotations

import pytest

import golden
from conftest import BLAST_DIR, METADATA_DIR, N_GENOMES

pytestmark = [pytest.mark.genomes, pytest.mark.tools]


def _hit_pairs(path):
    pairs = set()
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 12:
                pairs.add((fields[0], fields[1]))
    return pairs


@pytest.mark.parametrize("name", ["thereand.blast", "backagain.blast"])
def test_start_reproduces_the_committed_blast(pipeline, name):
    """The BLAST fixture the fast tier replays must match a real run.

    Compared on query/subject pairs rather than whole lines: bitscores and e-values move
    between BLAST versions, but which sequences hit which should not.
    """
    fresh = _hit_pairs(pipeline.blast / name)
    committed = _hit_pairs(BLAST_DIR / name)
    assert fresh == committed, (
        f"{name} drifted from tests/data/blast/{name}\n"
        f"  only in fresh run: {sorted(fresh - committed)[:5]}\n"
        f"  only in fixture  : {sorted(committed - fresh)[:5]}\n"
        "If this is expected (new BLAST version, changed inputs), regenerate the "
        "fixtures and the golden files together.")


def test_committed_metadata_covers_every_blast_protein(pipeline):
    """The trimmed metadata must still describe every protein the BLAST mentions."""
    referenced = set()
    for name, column in [("thereand.blast", 1), ("backagain.blast", 0)]:
        for query, subject in _hit_pairs(BLAST_DIR / name):
            value = subject if column == 1 else query
            if "__" in value:
                referenced.add(value)

    for name in ("genome_names.tsv", "genome_functions.tsv"):
        with open(METADATA_DIR / name) as fh:
            next(fh)
            have = {line.split("\t", 1)[0] for line in fh}
        assert not (referenced - have), f"{name} is missing {sorted(referenced - have)[:5]}"


def test_full_run_matches_the_analyze_snapshots(request, pipeline):
    """Running from genomes must give the same analysis as replaying the BLAST."""
    golden.assert_matches(
        request, "copy_count_matrix.tsv",
        (pipeline.analysis / "copy_count_matrix.tsv").read_text())
    golden.assert_matches(
        request, "expansion_summary.tsv",
        (pipeline.analysis / "expansion_summary.tsv").read_text())


def test_genome_count_survives_the_whole_pipeline(pipeline):
    rows = (pipeline.analysis / "copy_count_matrix.tsv").read_text().splitlines()
    assert len(rows) - 1 == N_GENOMES
