"""`trees` with MIBiG references — the EvoMining prediction step.

Needs the genomes, the toolchain, and the MIBiG protein FASTA (downloaded separately;
see the README).

MIBiG is what makes this tier worth having. `find_evomining_predictions` walks outward
from each MIBiG reference leaf until it hits a BBH leaf, marking the normal leaves it
passes as predictions — so with no MIBiG there are no seeds, no walk, and the green
"EvoMining prediction" class is unreachable. That is the actual scientific output of the
pipeline, so it is tested here rather than left to a no-MIBiG smoke run.
"""

from __future__ import annotations

import pytest

from evomining import trees

import golden

pytestmark = [pytest.mark.genomes, pytest.mark.tools, pytest.mark.mibig]

FAMILY_DIR = "3PGA_AMINOACIDS_1"

GREEN = trees.CLASSIFICATION_LABELS["evomining_prediction"]
MIBIG_REF = trees.CLASSIFICATION_LABELS["mibig_ref"]
CONSERVED = trees.CLASSIFICATION_LABELS["conserved"]
SEED = trees.CLASSIFICATION_LABELS["central"]


def _annotations(mibig_run):
    rows = []
    lines = (mibig_run.trees / FAMILY_DIR / "annotations.tsv").read_text().splitlines()
    assert lines[0].split("\t") == ["leaf_id", "classification", "color"]
    for line in lines[1:]:
        if line:
            rows.append(line.split("\t"))
    return rows


def test_mibig_blast_is_produced(mibig_run):
    blast = mibig_run.analysis / "expanded_vs_mibig.blast"
    assert blast.exists() and blast.stat().st_size > 0
    assert (mibig_run.analysis / "expanded_queries.faa").exists()


def test_classified_proteins_gain_mibig_columns(mibig_run):
    """With --mibig, the classification table must actually carry hits.

    Without this, a silently failing MIBiG BLAST would still produce a well-formed table
    full of NA and every other assertion would pass.
    """
    lines = (mibig_run.analysis / "classified_proteins.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    idx = header.index("mibig_hit")
    hits = [l.split("\t")[idx] for l in lines[1:] if l]
    assert any(h != "NA" for h in hits), "no protein got a MIBiG hit"


def test_tree_contains_mibig_reference_leaves(mibig_run):
    labels = [label for _, label, _ in _annotations(mibig_run)]
    assert labels.count(MIBIG_REF) > 0, "no MIBiG reference leaves were placed"
    assert labels.count(SEED) == 4, "the 4 Central DB seeds should each be a leaf"


def test_evomining_predictions_are_made(mibig_run):
    """The green class must be populated, and only from genome leaves.

    A prediction is a genome copy that sits between a MIBiG reference and the nearest
    conserved (BBH) leaf. Seeds and MIBiG references can never be predictions, so their
    appearance in the green set would mean the walk is mislabelling leaves.
    """
    rows = _annotations(mibig_run)
    green = [leaf for leaf, label, _ in rows if label == GREEN]
    assert green, "no EvoMining predictions were made"

    for leaf in green:
        assert not leaf.startswith("MIBIG|"), f"MIBiG reference marked green: {leaf}"
        assert not leaf.startswith("CENTRAL|"), f"Central seed marked green: {leaf}"


def test_prediction_and_conserved_classes_are_disjoint(mibig_run):
    """A leaf is one thing. The walk stops at conserved leaves, so it cannot mark one."""
    by_leaf = {leaf: label for leaf, label, _ in _annotations(mibig_run)}
    green = {l for l, lab in by_leaf.items() if lab == GREEN}
    conserved = {l for l, lab in by_leaf.items() if lab == CONSERVED}
    assert not (green & conserved)


def test_tree_leaves_match_the_sequences_that_went_in(mibig_run):
    """Every aligned sequence appears as a leaf, and nothing else does.

    Catches a leaf-id sanitisation change silently dropping or renaming members.
    """
    fam = mibig_run.trees / FAMILY_DIR
    sequences = {line[1:].strip()
                 for line in (fam / "sequences.faa").read_text().splitlines()
                 if line.startswith(">")}
    annotated = {leaf for leaf, _, _ in _annotations(mibig_run)}
    assert annotated == sequences


def test_classification_counts(request, mibig_run):
    """Snapshot how many leaves land in each class — but not the tree topology.

    Topology depends on MUSCLE and FastTree versions rather than on anything this
    codebase decides, so pinning it would produce failures that say nothing about
    evomining-py. Which leaf gets which colour is decided here, so that is pinned.

    Note the leaf *order* in sequences.faa is not reproducible: `collect_family_sequences`
    iterates a set of MIBiG subjects, so ordering varies with Python's per-process hash
    seed. The resulting classifications were stable across repeated runs on this dataset,
    but counts are snapshotted rather than the file itself for that reason.
    """
    counts = {}
    for _, label, colour in _annotations(mibig_run):
        counts[(label, colour)] = counts.get((label, colour), 0) + 1

    report = "classification\tcolor\tn_leaves\n" + "".join(
        f"{label}\t{colour}\t{n}\n" for (label, colour), n in sorted(counts.items()))
    golden.assert_matches(request, f"tree_{FAMILY_DIR}_classifications.tsv", report)


def test_mibig_leaves_are_added_in_sorted_order(mibig_run):
    """Guards the fix for non-reproducible leaf ordering.

    `collect_family_sequences` used to iterate a raw set of MIBiG subjects, so the record
    order in sequences.faa varied with Python's per-process hash seed -- meaning two runs
    over identical data handed MUSCLE different input and could produce different tree
    topologies. Sorting fixed it; this checks the ordering it guarantees.
    """
    fam = mibig_run.trees / FAMILY_DIR
    leaves = [line[1:].strip()
              for line in (fam / "sequences.faa").read_text().splitlines()
              if line.startswith(">")]
    mibig = [leaf for leaf in leaves if leaf.startswith("MIBIG")]
    assert mibig, "no MIBiG references in the alignment input"
    assert mibig == sorted(mibig), "MIBiG leaves are not in sorted order"
    assert leaves[-len(mibig):] == mibig, "MIBiG references should be appended last"


def test_colours_and_labels_come_from_the_known_sets(mibig_run):
    valid_colours = set(trees.CLASSIFICATION_COLORS.values())
    valid_labels = set(trees.CLASSIFICATION_LABELS.values())
    for leaf, label, colour in _annotations(mibig_run):
        assert colour in valid_colours, f"{leaf}: unexpected colour {colour}"
        assert label in valid_labels, f"{leaf}: unexpected label {label}"


def test_only_expanded_families_get_trees(mibig_run):
    built = {p.name for p in mibig_run.trees.iterdir() if p.is_dir()}
    assert built == {FAMILY_DIR}
    assert "# Total skipped: 0" in (mibig_run.trees / "skipped_families.txt").read_text()


def test_optional_exports_are_written(mibig_run):
    fam = mibig_run.trees / FAMILY_DIR
    itol = (fam / "itol_colors.txt").read_text().splitlines()
    assert itol[:3] == ["TREE_COLORS", "SEPARATOR TAB", "DATA"]

    micro = (fam / "microreact.csv").read_text().splitlines()
    assert micro[0] == "id,classification,classification__colour"
    # iTOL writes two rows per leaf (branch + label), Microreact one.
    assert len(micro) - 1 == len(itol[3:]) // 2


def test_tree_and_svg_exist(mibig_run):
    fam = mibig_run.trees / FAMILY_DIR
    tree = fam / "tree_rooted.nwk"
    if not tree.exists():
        tree = fam / "tree.nwk"
    assert tree.exists() and tree.stat().st_size > 0
    assert (fam / "tree.svg").stat().st_size > 0
