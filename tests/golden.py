"""Golden-file (snapshot) comparison.

The suite is a characterization suite: evomining-py's current output on the los17 example
genomes *is* the reference.  These helpers compare a fresh run against the committed
snapshot and, on request, rewrite it.

The point is not that the snapshot is correct -- it is that any change to it is visible
and deliberate.  When a diff appears, the question to ask is "did I mean to change this?",
and if the answer is yes::

    pytest --update-golden

then review the resulting diff in git before committing it.  A snapshot updated without
reading the diff is worse than no snapshot at all, because it launders a regression into
the baseline.
"""

from __future__ import annotations

import difflib
from pathlib import Path

import pytest

GOLDEN = Path(__file__).resolve().parent / "golden"

#: Lines matching these are dropped before comparison -- they vary run to run without
#: meaning anything changed.
_VOLATILE_PREFIXES = ("# generated", "# created")


def normalise(text):
    """Strip trailing whitespace and volatile lines so diffs show real changes only."""
    lines = []
    for line in text.splitlines():
        if line.startswith(_VOLATILE_PREFIXES):
            continue
        lines.append(line.rstrip())
    return "\n".join(lines).rstrip() + "\n"


def assert_matches(request, name, actual):
    """Compare `actual` against the snapshot `name`, or rewrite it under --update-golden."""
    path = GOLDEN / name
    actual = normalise(actual)

    if request.config.getoption("--update-golden"):
        path.parent.mkdir(parents=True, exist_ok=True)
        changed = not path.exists() or normalise(path.read_text()) != actual
        path.write_text(actual)
        if changed:
            print(f"\n[golden] updated {name}")
        return

    if not path.exists():
        pytest.fail(
            f"no snapshot at {path.relative_to(GOLDEN.parent)}.\n"
            f"Create it with:  pytest --update-golden")

    expected = normalise(path.read_text())
    if expected == actual:
        return

    diff = "\n".join(difflib.unified_diff(
        expected.splitlines(), actual.splitlines(),
        fromfile=f"golden/{name}", tofile="this run", lineterm="", n=2))
    pytest.fail(
        f"output changed against golden/{name}\n\n{diff}\n\n"
        f"If this change is intended, re-run with --update-golden and review the git diff.")


def tsv_without_columns(path, drop):
    """Read a TSV back as text with `drop` columns removed.

    Used to snapshot tables while excluding fields that legitimately move between runs
    or tool versions (BLAST bitscores and e-values, for instance) without giving up on
    comparing everything else in the row.
    """
    lines = Path(path).read_text().splitlines()
    header = lines[0].split("\t")
    keep = [i for i, col in enumerate(header) if col not in drop]
    out = []
    for line in lines:
        fields = line.split("\t")
        out.append("\t".join(fields[i] for i in keep if i < len(fields)))
    return "\n".join(out) + "\n"
