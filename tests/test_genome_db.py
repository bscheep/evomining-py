"""`generate-genome-db` on the los17 GenBank set.

Needs ``tests/data/los17_genomes/`` (233 MB, not committed -- see ``accessions.txt``).
No external tools: this stage is pure Biopython.
"""

from __future__ import annotations

import csv
from types import SimpleNamespace

import pytest

from evomining import genomedb
from evomining.io.genbank import load_genbank

import golden
from conftest import N_GENOMES

pytestmark = pytest.mark.genomes


@pytest.fixture(scope="session")
def genome_db(genomes_dir, tmp_path_factory):
    out = tmp_path_factory.mktemp("genomedb") / "evomining_db"
    genomedb.run(SimpleNamespace(
        input_dir=str(genomes_dir), output_dir=str(out),
        lists=None, list_dir=None, names=None, keep_pseudogenes=False))
    return out


def _proteins(genome_db):
    return [line[1:].strip()
            for line in (genome_db / "GENOMES.fasta").read_text().splitlines()
            if line.startswith(">")]


def test_per_genome_protein_counts(request, genome_db):
    """Snapshot how many proteins each genome contributes.

    This is the number every downstream count depends on, and it is the first thing a
    change in the parser (pseudogene handling, isoform collapsing, id fallback) would
    move.
    """
    counts = {}
    with open(genome_db / "genome_names.tsv") as fh:
        next(fh)
        for line in fh:
            stem = line.split("\t", 1)[0].split("__", 1)[0]
            counts[stem] = counts.get(stem, 0) + 1

    report = "genome_stem\tproteins\n" + "".join(
        f"{stem}\t{counts[stem]}\n" for stem in sorted(counts))
    golden.assert_matches(request, "genome_db_protein_counts.tsv", report)


def test_genome_set_is_complete(genome_db):
    stems = {p.split("__", 1)[0] for p in _proteins(genome_db)}
    assert len(stems) == N_GENOMES
    assert all(s.startswith("GCF_") for s in stems), sorted(stems)[:3]


def test_composite_ids_are_unique_and_well_formed(genome_db):
    """`<stem>__<locus_tag>` is the namespace every later artifact is keyed on.

    A duplicate would silently merge two proteins downstream, so uniqueness is checked
    rather than assumed.
    """
    proteins = _proteins(genome_db)
    assert len(proteins) == len(set(proteins)), "composite protein ids must be unique"
    for pid in proteins:
        stem, sep, locus = pid.partition("__")
        assert sep and stem and locus, pid
        assert " " not in pid, "BLAST truncates ids at whitespace"


def test_artifacts_line_up_with_each_other(genome_db):
    """The three outputs must describe exactly the same protein set."""
    proteins = set(_proteins(genome_db))
    for name, header in [("genome_names.tsv", "protein_id\tgenome_name"),
                         ("genome_functions.tsv", "protein_id\tfunction")]:
        lines = (genome_db / name).read_text().splitlines()
        assert lines[0] == header
        ids = {line.split("\t", 1)[0] for line in lines[1:]}
        assert ids == proteins, f"{name} does not match GENOMES.fasta"


def test_fasta_is_wrapped_at_60_columns(genome_db):
    for line in (genome_db / "GENOMES.fasta").read_text().splitlines():
        if not line.startswith(">"):
            assert len(line) <= 60


def test_every_protein_has_a_name_and_function(genome_db):
    with open(genome_db / "genome_functions.tsv") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))[1:]
    assert all(len(r) == 2 and r[1] for r in rows)

    hypothetical = sum(1 for r in rows if r[1] == "hypothetical protein")
    assert hypothetical < len(rows) * 0.5, (
        f"{hypothetical}/{len(rows)} products are 'hypothetical protein'")


def test_pseudogenes_are_dropped_by_default(genomes_dir):
    """Pseudogenes are not functional enzymes and would inflate copy counts.

    Checked on one genome rather than the whole set, because parsing 233 MB twice is
    slow and the flag is parser-wide.
    """
    sample = sorted(genomes_dir.glob("*.gbff"))[0]
    dropped = load_genbank(sample)
    kept = load_genbank(sample, keep_pseudo=True)
    assert kept.n_genes >= dropped.n_genes
    assert dropped.n_genes > 0
