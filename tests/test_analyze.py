"""`analyze` on the los17 example: snapshots plus the invariants worth stating outright.

Runs from committed BLAST output, so this whole module needs no genomes and no external
tools.  It is the fast tier and the one that will catch most regressions, because every
number the pipeline reports downstream is derived here.
"""

from __future__ import annotations

import csv

import pytest

from evomining import analysis

import golden

PGDH = "3PGA_AMINOACIDS|1"   # phosphoglycerate dehydrogenase -- the expanded family
PSAT = "3PGA_AMINOACIDS|2"   # phosphoserine aminotransferase -- the flat control


def _read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------------------
#  Snapshots
# ---------------------------------------------------------------------------

def test_copy_count_matrix(request, analysis_out):
    golden.assert_matches(
        request, "copy_count_matrix.tsv",
        (analysis_out / "copy_count_matrix.tsv").read_text())


def test_expansion_summary(request, analysis_out):
    golden.assert_matches(
        request, "expansion_summary.tsv",
        (analysis_out / "expansion_summary.tsv").read_text())


def test_classified_proteins(request, analysis_out):
    """Snapshot the classification table without the BLAST-derived numeric columns.

    `evalue` and `bitscore` come straight from BLAST and shift with its version, so
    pinning them would make the snapshot fail for reasons that have nothing to do with
    evomining-py. Everything that evomining-py actually decides -- family, genome,
    protein, function, BBH, classification, copy number, expansion -- is compared.
    """
    golden.assert_matches(
        request, "classified_proteins.tsv",
        golden.tsv_without_columns(
            analysis_out / "classified_proteins.tsv",
            drop={"evalue", "bitscore", "mibig_pident", "mibig_evalue"}))


# ---------------------------------------------------------------------------
#  Invariants -- stated explicitly so a snapshot update cannot quietly break them
# ---------------------------------------------------------------------------

def test_every_genome_appears_exactly_once(analysis_out, genome_names):
    rows = _read(analysis_out / "copy_count_matrix.tsv")
    stems = [r["genome_stem"] for r in rows]
    assert len(stems) == len(set(stems)) == len(genome_names)


def test_the_two_families_behave_as_the_example_intends(analysis_out):
    """PGDH varies and expands; PSAT is single-copy everywhere.

    This is the biological point of the example dataset, and it holds regardless of the
    exact expansion threshold -- so it is asserted separately from the snapshot.
    """
    rows = {r["enzyme_family"]: r for r in _read(analysis_out / "expansion_summary.tsv")}
    assert set(rows) == {PGDH, PSAT}

    assert int(rows[PGDH]["n_expanded_genomes"]) > 0, "PGDH must be expanded somewhere"
    assert int(rows[PSAT]["n_expanded_genomes"]) == 0, "PSAT is the flat control"
    assert float(rows[PSAT]["std"]) == 0.0, "PSAT is single-copy in every genome"
    assert int(rows[PGDH]["max_copies"]) > int(rows[PSAT]["max_copies"])


def test_expanded_genomes_are_above_their_family_threshold(analysis_out):
    """The expansion column must agree with the counts it claims to summarise."""
    counts = {}
    for row in _read(analysis_out / "copy_count_matrix.tsv"):
        for family in (PGDH, PSAT):
            counts[(family, row["genome_name"])] = int(row[family])

    for row in _read(analysis_out / "expansion_summary.tsv"):
        family, threshold = row["enzyme_family"], float(row["threshold"])
        listed = [g for g in row["expanded_genomes"].split("; ") if g]
        assert len(listed) == int(row["n_expanded_genomes"])
        for entry in listed:
            name, _, tail = entry.rpartition("(")
            count = int(tail.rstrip(")"))
            assert count > threshold, f"{family} {name}: {count} is not > {threshold}"
            assert counts[(family, name.strip())] == count


def test_classification_is_consistent_with_bbh(analysis_out):
    """With no antiSMASH mapping, every copy is either `conserved` (BBH) or `normal`.

    `classify_proteins` folds antiSMASH membership into the label, so this pins the
    no-antiSMASH case and would catch a stray cyan/purple label appearing from nowhere.
    """
    for row in _read(analysis_out / "classified_proteins.tsv"):
        assert row["in_antismash"] == "False"
        expected = "conserved" if row["is_bbh"] == "True" else "normal"
        assert row["classification"] == expected, row


def test_bbh_proteins_are_a_subset_of_counted_proteins(analysis_out, forward_hits,
                                                       reverse_hits, genome_names):
    """A protein can only be a conserved copy if it was counted as a copy at all."""
    _, members = analysis.compute_copy_counts(forward_hits, genome_names, method="multi")
    counted = {m["protein_id"]
               for genomes in members.values() for ms in genomes.values() for m in ms}
    bbh = analysis.compute_bbh(forward_hits, reverse_hits)

    flagged = {row["protein_id"] for row in _read(analysis_out / "classified_proteins.tsv")
               if row["is_bbh"] == "True"}
    assert flagged <= counted
    assert flagged <= bbh


# ---------------------------------------------------------------------------
#  Copy-count method
# ---------------------------------------------------------------------------

def test_best_method_never_exceeds_multi(forward_hits, genome_names):
    """`--ccm best` assigns each protein to one family; `multi` counts it in all of them.

    So `best` totals can never exceed `multi` totals. This is the invariant that makes
    the two modes meaningfully different rather than accidentally equal.
    """
    multi, _ = analysis.compute_copy_counts(forward_hits, genome_names, method="multi")
    best, _ = analysis.compute_copy_counts(forward_hits, genome_names, method="best")

    total_multi = sum(c for fam in multi.values() for c in fam.values())
    total_best = sum(c for fam in best.values() for c in fam.values())
    assert total_best <= total_multi

    for family, genomes in best.items():
        for stem, count in genomes.items():
            assert count <= multi[family][stem]


def test_unknown_copy_count_method_is_rejected(forward_hits, genome_names):
    with pytest.raises(ValueError, match="Unknown copy count method"):
        analysis.compute_copy_counts(forward_hits, genome_names, method="bogus")
