#!/usr/bin/env python3
"""
evomining.trees

Build colored phylogenetic trees for expanded EvoMining enzyme families,
following the EvoMining methodology (Selem-Mojica et al. 2019):

  1. Collect sequences: Central DB seeds + genome homologs + MIBiG hits
  2. Align with MUSCLE (fallback MAFFT)
  3. Trim with trimAl -automated1
  4. Build tree with FastTree 2.1 (approx ML)
  5. Root on BBH (conserved metabolism) leaf
  6. Identify EvoMining predictions via nw_clade (same algorithm as color_tree.pl)
  7. Color and render an SVG; write annotations.tsv

Composite-ID version. Genome protein sequences are read from the single
GENOMES.fasta written by `generate-genome-db`, keyed on composite
<genome_stem>__<locus_tag> IDs. There is no GENOMES/ directory of per-genome
.faa, no fig|/RAST id reconstruction, and no Corason_Rast.IDs indirection.
classified_proteins.tsv supplies protein_id, genome_name, function, and the
classification used for leaf colouring.

Clade / iTOL metadata is intentionally NOT produced here: clade assignment is a
study-specific overlay, not part of EvoMining. trees emits the tree, its SVG,
and annotations.tsv; iTOL-ready metadata files are the job of a separate tool.

Usage:
  evomining trees --mibig /path/to/MiBIG_DB.faa
  evomining trees --classified evomining_analysis/classified_proteins.tsv \\
                  --central-db enzymes_db.faa --genomes-fasta evomining_db/GENOMES.fasta \\
                  --families symb2_families.txt -o evomining_trees
"""

import argparse
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from .external import Tools


# ---------------------------------------------------------------------------
#  Constants
# ---------------------------------------------------------------------------

CLASSIFICATION_COLORS = {
    "conserved":            "#FF0000",  # red        — BBH conserved metabolism
    "antismash":            "#00CED1",  # cyan       — antiSMASH predicted
    "transition":           "#800080",  # purple     — both BBH and antiSMASH
    "mibig_hit":            "#0000FF",  # blue       — genome copy with MIBiG BLAST hit
    "evomining_prediction": "#a5bb47",  # green      — EvoMining prediction (post-tree)
    "normal":               "#808080",  # grey       — all other copies
    "central":              "#FFA500",  # gold       — Central DB seed
    "mibig_ref":            "#4169E1",  # royal blue — MIBiG reference sequence
}

# Map internal classification names to original EvoMining labels
CLASSIFICATION_LABELS = {
    "central":             "Seed Enzyme",
    "conserved":           "Central metabolism",
    "antismash":           "Secondary metabolism (antiSMASH)",
    "transition":          "Transition Enzyme (antiSMASH and Central metabolism)",
    "mibig_ref":           "Recruited Enzyme (MIBiG)",
    "evomining_prediction":"Secondary metabolism (EvoMining Hit)",
    "normal":              "expansion",
    "mibig_hit":           "Recruited Enzyme (MIBiG)",
}


# ---------------------------------------------------------------------------
#  EvoMining prediction algorithm — mirrors color_tree.pl exactly (UNCHANGED)
# ---------------------------------------------------------------------------

def find_evomining_predictions(tree_file, tree_seqs, tools, max_depth=30):
    """
    Identify EvoMining predictions using the same algorithm as color_tree.pl.

    For each MIBiG reference leaf:
      1. nw_clade -c N expands the clade at increasing depth
      2. Skip non-genome leaves (central)
      3. Red (BBH) leaf -> clear candidates, mark seed exhausted, stop
      4. Normal leaves not red/cyan/blue -> green predictions

    Returns: set of tree_ids to reclassify as 'evomining_prediction'
    """
    mibig_leaves     = set()  # mibig_ref seeds
    bbh_leaves       = set()  # conserved + transition (red boundary)
    antismash_leaves = set()
    non_genome_leaves = set()  # central (skip in clade loop)

    for tree_id, (_, _, classification) in tree_seqs.items():
        if classification == "mibig_ref":
            mibig_leaves.add(tree_id)
        elif classification in ("conserved", "transition"):
            bbh_leaves.add(tree_id)
        elif classification == "antismash":
            antismash_leaves.add(tree_id)
        elif classification == "central":
            non_genome_leaves.add(tree_id)

    if not mibig_leaves:
        return set()

    predictions  = set()
    exhausted_kr = set()

    for depth in range(0, max_depth + 1):
        for kr_leaf in list(mibig_leaves - exhausted_kr):
            r = run_cmd(tools, "nw_clade",
                        ["-c", str(depth), str(tree_file), kr_leaf], quiet=True)
            if r.returncode != 0 or not r.stdout.strip():
                continue

            r2 = run_cmd(tools, "nw_labels", ["-I", "-"], stdin=r.stdout, quiet=True)
            if r2.returncode != 0:
                continue

            clade_leaves = [l.strip() for l in r2.stdout.strip().split("\n") if l.strip()]

            cadena_kr   = []
            flag_cadena = False
            hit_red     = False

            for leaf in clade_leaves:
                # Skip non-genome leaves (mirrors: if($x !~ /\|/){ next; })
                if leaf in non_genome_leaves:
                    continue
                # Red boundary (mirrors: $cadenaKR=''; delete $hashKR{$y}; last)
                if leaf in bbh_leaves:
                    cadena_kr   = []
                    flag_cadena = False
                    hit_red     = True
                    break
                # Normal, not cyan, not blue -> candidate
                if leaf not in antismash_leaves and leaf not in mibig_leaves:
                    cls = tree_seqs.get(leaf, (None, None, "unknown"))[2]
                    if cls == "normal":
                        cadena_kr.append(leaf)
                        flag_cadena = True

            if hit_red:
                exhausted_kr.add(kr_leaf)
                continue

            # mirrors: if($flagcadenaKR==1){ print KR/KRR }
            if flag_cadena:
                for leaf in cadena_kr:
                    predictions.add(leaf)

    return predictions


# ---------------------------------------------------------------------------
#  SVG colorizer
# ---------------------------------------------------------------------------

def colorize_svg(svg_content, tree_seqs, rename_file):
    """Post-process nw_display SVG to color leaf labels by classification.

    nw_display renders leaf labels with underscores converted to spaces, so the
    label in the SVG text element does not match the underscored rename_map label
    verbatim. We therefore match on the space-substituted form of each label.
    """
    rename_map = {}
    with open(rename_file) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                rename_map[parts[0]] = parts[1]

    label_to_color = {}
    for tree_id, (_, _, classification) in tree_seqs.items():
        label = rename_map.get(tree_id, tree_id)
        color = CLASSIFICATION_COLORS.get(classification, "#808080")
        # nw_display shows underscores as spaces in the rendered SVG text.
        svg_label = label.replace("_", " ")
        label_to_color[svg_label] = color

    # Longest labels first, so a short label isn't matched as a prefix of a long one.
    for label, color in sorted(label_to_color.items(), key=lambda kv: -len(kv[0])):
        escaped = re.escape(label)
        # nw_display emits <text ...>LABEL</text>; wrap the label text in a
        # colored tspan. Match the label as the full text content of the element.
        svg_content = re.sub(
            f'(<text[^>]*>){escaped}(</text>)',
            f'\\1<tspan fill="{color}">{label}</tspan>\\2',
            svg_content,
        )
    return svg_content


# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------

def run_cmd(tools, tool, args, stdin=None, quiet=False):
    """Run an external tool via the resolver (PATH, or a configured env)."""
    return tools.run(tool, args, stdin=stdin, quiet=quiet)


def read_fasta(path):
    """Read FASTA -> dict: header_id -> (full_header, sequence)."""
    seqs = {}
    current_id = current_header = None
    current_seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = (current_header, "".join(current_seq))
                current_header = line
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        seqs[current_id] = (current_header, "".join(current_seq))
    return seqs


def write_fasta(seqs, path):
    """Write dict of id -> (header, seq, ...) to FASTA."""
    with open(path, "w") as fh:
        for sid, vals in seqs.items():
            fh.write(f"{vals[0]}\n")
            for i in range(0, len(vals[1]), 60):
                fh.write(vals[1][i:i+60] + "\n")


def load_genome_subset(genomes_fasta, wanted_ids):
    """Read GENOMES.fasta, returning {protein_id: (header, seq)} for wanted_ids only.

    protein_id is the leading whitespace-free token of each header, i.e. the
    composite <genome_stem>__<locus_tag>. Only the requested proteins are kept,
    so memory scales with the tree members, not the whole 2M-protein DB.
    """
    wanted = set(wanted_ids)
    seqs = {}
    cur_id = cur_header = None
    cur_seq = []
    with open(genomes_fasta) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_id in wanted and cur_seq:
                    seqs[cur_id] = (cur_header, "".join(cur_seq))
                cur_header = line
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id in wanted and cur_seq:
            seqs[cur_id] = (cur_header, "".join(cur_seq))
    return seqs


# ---------------------------------------------------------------------------
#  Load data
# ---------------------------------------------------------------------------

def load_classified(classified_file):
    """Load classified_proteins.tsv -> dict: family -> list of protein records."""
    families = defaultdict(list)
    with open(classified_file) as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            rec = dict(zip(header, parts))
            families[rec["family"]].append(rec)
    return dict(families)


# ---------------------------------------------------------------------------
#  Collect sequences for one family
# ---------------------------------------------------------------------------

def collect_family_sequences(family, members, central_seqs, genome_seqs, mibig_seqs, mibig_hits):
    """
    Collect all sequences for a family tree.

    Leaf ID format:
      Central DB:  CENTRAL|pathway|num|enzyme_name|source
      Genome:      {protein_id}|{genome_name}|{function}
      MIBiG ref:   MIBIG|{subject_id}

    Classification lives only in the tree_seqs tuple, never in the leaf ID.
    Genome sequences come from GENOMES.fasta (genome_seqs), keyed by the
    composite protein_id -- no per-genome .faa, no fig| reconstruction. The
    function is taken from classified_proteins.tsv (already backfilled from
    genome_functions.tsv by `analyze`).

    Returns dict: tree_id -> (header, sequence, classification)
    """
    tree_seqs = {}

    def sanitize_tree_id(tid):
        for ch in [",", "(", ")", ":", ";", "'", '"', "[", "]", " "]:
            tid = tid.replace(ch, "_")
        while "__" in tid:
            tid = tid.replace("__", "_")
        return tid.rstrip("_")

    # 1. Central DB seeds (matched by SUBSYSTEM|family_number)
    fam_parts = family.split("|")
    pathway   = fam_parts[0] if fam_parts else family
    fam_num   = fam_parts[1] if len(fam_parts) > 1 else ""

    for cid, (cheader, cseq) in central_seqs.items():
        cparts = cid.split("|")
        if len(cparts) >= 2 and cparts[0] == pathway and cparts[1] == fam_num:
            tree_id = sanitize_tree_id(f"CENTRAL|{cid}")
            tree_seqs[tree_id] = (f">{tree_id}", cseq, "central")

    # 2. Genome homologs (sequences from GENOMES.fasta by composite protein_id)
    for member in members:
        protein_id  = member["protein_id"]
        pre_class   = member["classification"]
        genome_name = member.get("genome_name", ".") or "."

        entry = genome_seqs.get(protein_id)
        if not entry:
            continue
        _, seq = entry

        func = (member.get("function", "unknown") or "unknown")[:60]
        short_name = genome_name[:60]

        # Leaf ID: {protein_id}|{genome_name}|{function} — no classification prefix
        tree_id = sanitize_tree_id(f"{protein_id}|{short_name}|{func}")
        tree_seqs[tree_id] = (f">{tree_id}", seq, pre_class)

    # 3. MIBiG reference sequences (all hits per genome protein, subject-level unique; bitscore >= 100)
    MAX_SEQ_LEN = 15000
    subjects_to_add = set()
    for m in members:
        subs = mibig_hits.get(m["protein_id"])
        if subs:
            subjects_to_add |= subs
    for mibig_subject in sorted(subjects_to_add): # sort for reproducibility between runs
        if mibig_subject in mibig_seqs:
            _, mseq = mibig_seqs[mibig_subject]
            if len(mseq) > MAX_SEQ_LEN:
                print(f"  WARNING: skipping MIBiG {mibig_subject} ({len(mseq)} aa, likely corrupted)")
                continue
            tree_id = sanitize_tree_id(f"MIBIG|{mibig_subject}")
            tree_seqs[tree_id] = (f">{tree_id}", mseq, "mibig_ref")

    return tree_seqs


# ---------------------------------------------------------------------------
#  Build tree for one family
# ---------------------------------------------------------------------------

def build_family_tree(family, tree_seqs, outdir, tools, threads=32):
    """
    Build phylogenetic tree for one enzyme family:
      1. Align with MUSCLE (fallback: MAFFT)
      2. Trim with trimAl -automated1
      3. Build tree with FastTree (approx ML)
      4. Root on BBH leaf (or Central DB seed if no BBH)

    The colored SVG is NOT rendered here; render_tree_svg() does that from run()
    AFTER find_evomining_predictions() reclassifies leaves, so green predictions
    are coloured correctly. annotations.tsv is likewise written later.

    Returns path to newick tree file (str), or None on failure.
    """
    safe_name = family.replace("|", "_")
    fam_dir   = Path(outdir) / safe_name
    fam_dir.mkdir(parents=True, exist_ok=True)

    input_fasta = fam_dir / "sequences.faa"
    write_fasta(tree_seqs, str(input_fasta))

    if len(tree_seqs) < 4:
        print(f"  {family}: only {len(tree_seqs)} sequences, skipping tree")
        return None

    # 1. Align
    aligned_fasta = fam_dir / "aligned.faa"
    r = run_cmd(tools, "muscle",
                ["-align", input_fasta, "-output", aligned_fasta, "-threads", threads],
                quiet=True)
    if not aligned_fasta.exists() or aligned_fasta.stat().st_size == 0:
        # MUSCLE v3 argument style, in case an older muscle is on PATH
        r = run_cmd(tools, "muscle", ["-in", input_fasta, "-out", aligned_fasta], quiet=True)
    if not aligned_fasta.exists() or aligned_fasta.stat().st_size == 0:
        print("  MUSCLE failed, trying MAFFT...")
        r = run_cmd(tools, "mafft", ["--auto", "--thread", threads, input_fasta])
        if r.returncode == 0 and r.stdout:
            aligned_fasta.write_text(r.stdout)
        else:
            print(f"  {family}: alignment failed (MUSCLE and MAFFT)")
            return None

    # 2. Trim
    trimmed_fasta = fam_dir / "trimmed.faa"
    r = run_cmd(tools, "trimal",
                ["-in", aligned_fasta, "-out", trimmed_fasta, "-automated1"])
    if not trimmed_fasta.exists() or trimmed_fasta.stat().st_size == 0:
        print(f"  trimAl failed, using untrimmed alignment")
        trimmed_fasta = aligned_fasta
    else:
        print(f"  trimAl: trimmed alignment written")

    # 3. Build tree
    tree_file = fam_dir / "tree.nwk"
    for binary in ("FastTree", "fasttree"):
        if not tools.available(binary):
            continue
        r = run_cmd(tools, binary, ["-quiet", trimmed_fasta])
        if r.returncode == 0 and r.stdout:
            tree_file.write_text(r.stdout)
            break

    if not tree_file.exists() or tree_file.stat().st_size == 0:
        print(f"  {family}: tree building failed")
        return None

    # 4. Root
    bbh_leaf = next(
        (tid for tid, (_, _, cls) in tree_seqs.items() if cls == "conserved"), None
    )
    root_leaf = bbh_leaf or next(
        (tid for tid, (_, _, cls) in tree_seqs.items() if cls == "central"), None
    )

    if root_leaf:
        rooted_tree = fam_dir / "tree_rooted.nwk"
        r = run_cmd(tools, "nw_reroot", [tree_file, root_leaf])
        if r.returncode == 0 and r.stdout.strip():
            rooted_tree.write_text(r.stdout)
            tree_file = rooted_tree
            root_type = "BBH" if bbh_leaf else "Central DB seed"
            print(f"  Rooted on {root_type}: {root_leaf[:60]}...")
        else:
            print(f"  WARNING: nw_reroot failed, using unrooted tree")
    else:
        print(f"  No BBH or Central seed found for rooting, using unrooted tree")

    # NOTE: SVG rendering is deliberately NOT done here. It is done by
    # render_tree_svg() from run(), AFTER find_evomining_predictions() has
    # reclassified leaves, so green predictions are coloured correctly.
    return str(tree_file)


def render_tree_svg(family, tree_file, tree_seqs, outdir, tools):
    """
    Render the colored SVG for a family tree.

    Called from run() AFTER find_evomining_predictions() has reclassified
    leaves, so 'evomining_prediction' leaves render green. Writes rename_map.tsv,
    tree_renamed.nwk, and tree.svg into the family directory.
    """
    fam_dir = Path(outdir) / family.replace("|", "_")

    # Genome leaf ID: parts[0]=protein_id, parts[1]=genome_name, parts[2]=func
    rename_file = fam_dir / "rename_map.tsv"
    with open(rename_file, "w") as fh:
        for tree_id, (_, _, classification) in tree_seqs.items():
            parts = tree_id.split("|")
            if classification == "central":
                label = f"SEED_{parts[3]}_{parts[4]}" if len(parts) >= 5 else (
                    f"SEED_{parts[3]}" if len(parts) >= 4 else
                    f"SEED_{parts[-1]}" if len(parts) >= 2 else tree_id
                )
            elif classification == "mibig_ref":
                label = f"{parts[6]}_{parts[1]}" if len(parts) >= 7 else (
                    parts[1] if len(parts) >= 2 else tree_id
                )
            else:
                # Genome: parts[0]=protein_id, parts[1]=genome_name, parts[2]=func
                label = f"{parts[1]}_{parts[2]}" if len(parts) >= 3 else (
                    parts[1] if len(parts) >= 2 else tree_id
                )
            fh.write(f"{tree_id}\t{label}\n")

    renamed_tree = fam_dir / "tree_renamed.nwk"
    r = run_cmd(tools, "nw_rename", [tree_file, rename_file])
    if r.returncode == 0 and r.stdout:
        renamed_tree.write_text(r.stdout)
        r = run_cmd(tools, "nw_display", ["-s", "-w", "1200", renamed_tree])
        if r.returncode == 0 and r.stdout:
            svg_content = colorize_svg(r.stdout, tree_seqs, rename_file)
            (fam_dir / "tree.svg").write_text(svg_content)


# ---------------------------------------------------------------------------
#  Annotations
# ---------------------------------------------------------------------------

def write_annotations(tree_seqs, outdir, family):
    """
    Write annotations.tsv (leaf_id, classification label, color).
    Must be called AFTER find_evomining_predictions() has reclassified leaves
    so that 'evomining_prediction' entries are correctly captured.
    """
    fam_dir         = Path(outdir) / family.replace("|", "_")
    annotation_file = fam_dir / "annotations.tsv"
    with open(annotation_file, "w") as fh:
        fh.write("leaf_id\tclassification\tcolor\n")
        for tree_id, (_, _, classification) in tree_seqs.items():
            color = CLASSIFICATION_COLORS.get(classification, "#808080")
            label = CLASSIFICATION_LABELS.get(classification, classification)
            fh.write(f"{tree_id}\t{label}\t{color}\n")


def write_itol_colorstrip(tree_seqs, outdir, family):
    """
    Write an iTOL TREE_COLORS file: colour each leaf's branch and label by its
    classification, using CLASSIFICATION_COLORS. Opt-in via `--itol`.

    Keyed on the (sanitized) leaf tree_id, so it lines up with tree.nwk and
    annotations.tsv. Must be called AFTER find_evomining_predictions() so
    'evomining_prediction' leaves come out green. Upload tree.nwk to iTOL, then
    drag this file onto the tree.
    """
    fam_dir   = Path(outdir) / family.replace("|", "_")
    itol_file = fam_dir / "itol_colors.txt"
    with open(itol_file, "w") as fh:
        fh.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")
        for tree_id, (_, _, classification) in tree_seqs.items():
            color = CLASSIFICATION_COLORS.get(classification, "#808080")
            fh.write(f"{tree_id}\tbranch\t{color}\tnormal\t2\n")
            fh.write(f"{tree_id}\tlabel\t{color}\tnormal\t1\n")


def write_microreact_csv(tree_seqs, outdir, family):
    """
    Write a Microreact metadata CSV: id, classification, classification__colour.

    Opt-in via `--microreact`. The `id` column is the raw (sanitized) leaf
    tree_id, so it matches the leaf names in tree_rooted.nwk / tree.nwk exactly.

    Microreact colours a column X by an adjacent column named 'X__colour', so the
    hex colours live in 'classification__colour' and Microreact applies them to
    the 'classification' values automatically on upload.

    Must be called AFTER find_evomining_predictions() so 'evomining_prediction'
    leaves come out green. Load into microreact.org together with the .nwk tree.
    """
    fam_dir = Path(outdir) / family.replace("|", "_")
    csv_file = fam_dir / "microreact.csv"
    with open(csv_file, "w") as fh:
        fh.write("id,classification,classification__colour\n")
        for tree_id, (_, _, classification) in tree_seqs.items():
            color = CLASSIFICATION_COLORS.get(classification, "#808080")
            label = CLASSIFICATION_LABELS.get(classification, classification)
            # Quote fields to be safe against commas in labels.
            fh.write(f'"{tree_id}","{label}","{color}"\n')


# ---------------------------------------------------------------------------
#  HTML index
# ---------------------------------------------------------------------------

def write_tree_index_html(families_built, families_skipped, outdir):
    """Write index.html linking to all tree output files."""
    index_file = Path(outdir) / "index.html"
    html = """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>EvoMining Phylogenetic Trees</title>
<style>
  body { font-family: 'Courier New', monospace; background: #1a1a2e; color: #e0e0e0; padding: 20px; }
  h1 { color: #e94560; } h2 { color: #e94560; margin-top: 30px; }
  .legend { display: flex; gap: 15px; margin: 15px 0; font-size: 0.85em; flex-wrap: wrap; }
  .legend-item { display: flex; align-items: center; gap: 5px; }
  .legend-box { width: 14px; height: 14px; border: 1px solid #444; }
  table { border-collapse: collapse; width: 100%; margin-top: 15px; }
  th { background: #16213e; color: #e94560; padding: 8px 12px; text-align: left; border: 1px solid #333; }
  td { padding: 6px 12px; border: 1px solid #222; }
  tr:hover { background: #16213e; }
  a { color: #7ec8e3; text-decoration: none; } a:hover { text-decoration: underline; }
  .stats { color: #888; font-size: 0.85em; } .skipped { color: #ff9966; }
</style>
</head>
<body>
<h1>EvoMining Phylogenetic Trees</h1>
""" + f'<p class="stats">{len(families_built)} enzyme families with expansion trees &mdash; {len(families_skipped)} skipped</p>' + """
<div class="legend">
  <div class="legend-item"><div class="legend-box" style="background:#FF0000"></div> Conserved (BBH)</div>
  <div class="legend-item"><div class="legend-box" style="background:#00CED1"></div> antiSMASH</div>
  <div class="legend-item"><div class="legend-box" style="background:#800080"></div> Transition</div>
  <div class="legend-item"><div class="legend-box" style="background:#0000FF"></div> MIBiG hit</div>
  <div class="legend-item"><div class="legend-box" style="background:#a5bb47"></div> EvoMining prediction</div>
  <div class="legend-item"><div class="legend-box" style="background:#FFA500"></div> Central DB seed</div>
  <div class="legend-item"><div class="legend-box" style="background:#808080"></div> Normal copy</div>
</div>
<table>
<tr><th>Enzyme Family</th><th>Sequences</th><th>Tree</th><th>SVG</th><th>Annotations</th></tr>
"""
    for fam, info in sorted(families_built.items()):
        safe = fam.replace("|", "_")
        # Link to the actual tree that was built — tree_rooted.nwk when rooting
        # succeeded, else tree.nwk. info["tree_file"] holds the real path.
        tree_name = Path(info["tree_file"]).name if info.get("tree_file") else "tree.nwk"
        html += (f'<tr><td>{fam} ({info["enzyme_name"]})</td><td>{info["n_seqs"]}</td>'
                 f'<td><a href="{safe}/{tree_name}">{tree_name}</a></td>'
                 f'<td><a href="{safe}/tree.svg">tree.svg</a></td>'
                 f'<td><a href="{safe}/annotations.tsv">annotations.tsv</a></td></tr>\n')

    if families_skipped:
        html += "</table>\n<h2>Skipped Families</h2>\n<table>\n<tr><th>Family</th><th>Sequences</th><th>Reason</th></tr>\n"
        for fam, info in sorted(families_skipped.items()):
            html += f'<tr class="skipped"><td>{fam}</td><td>{info["n_seqs"]}</td><td>{info["reason"]}</td></tr>\n'

    html += "</table>\n</body>\n</html>"
    index_file.write_text(html)


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def run(args, tools=None):
    """Build a colored phylogenetic tree per expanded enzyme family."""
    tools = tools or Tools()
    tools.require("trimal", "nw_clade", "nw_labels", "nw_reroot", "nw_rename", "nw_display")
    if not (tools.available("muscle") or tools.available("mafft")):
        tools.require("muscle")
    if not (tools.available("FastTree") or tools.available("fasttree")):
        tools.require("FastTree")

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    skip_log_file  = outdir / "skipped_families.txt"
    skip_log_lines = []

    def log_skip(family, n_seqs, reason):
        msg = f"SKIPPED\t{family}\t{n_seqs} seqs\t{reason}"
        print(f"  {msg}")
        skip_log_lines.append(msg)

    print("Loading classified proteins...")
    all_classified = load_classified(args.classified)
    print(f"  {len(all_classified)} enzyme families")

    print("Loading Central DB sequences...")
    central_seqs = read_fasta(args.central_db)
    print(f"  {len(central_seqs)} seed sequences")

    mibig_seqs = {}
    if args.mibig and Path(args.mibig).exists():
        print("Loading MIBiG sequences...")
        mibig_seqs = read_fasta(args.mibig)
        print(f"  {len(mibig_seqs)} MIBiG proteins")

    mibig_hits = {}
    if args.mibig_blast and Path(args.mibig_blast).exists():
        with open(args.mibig_blast) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 12:
                    qid, sid, bitscore = parts[0], parts[1], float(parts[11])
                    if bitscore < 100:   # $score2 floor; project copy-count threshold
                        continue
                    mibig_hits.setdefault(qid, set()).add(sid)
        n_subjects = len({s for subs in mibig_hits.values() for s in subs})
        print(f"  {len(mibig_hits)} queries with MIBiG hits, "
              f"{n_subjects} distinct subjects (all hits, bitscore >= 100)")

    # Determine families to process
    expanded_families = [
        family for family, members in all_classified.items()
        if any(str(m.get("is_expanded", "")).lower() == "true" for m in members)
    ]

    if args.families:
        if len(args.families) == 1 and args.families[0].endswith(".txt"):
            with open(args.families[0]) as fh:
                family_set = set(line.strip() for line in fh if line.strip())
        else:
            family_set = set(args.families)
        expanded_families = [f for f in expanded_families if f in family_set]

    print(f"\n{len(expanded_families)} families with expansions to process")

    # Load only the genome sequences these families need, once, from GENOMES.fasta.
    wanted_ids = {m["protein_id"] for fam in expanded_families for m in all_classified[fam]}
    print(f"Loading genome sequences for {len(wanted_ids)} proteins from GENOMES.fasta...")
    genome_seqs = load_genome_subset(args.genomes_fasta, wanted_ids)
    print(f"  {len(genome_seqs)} sequences loaded\n")

    families_built   = {}
    families_skipped = {}

    for i, family in enumerate(sorted(expanded_families)):
        members   = all_classified[family]
        safe_name = family.replace("|", "_")
        fam_dir   = Path(outdir) / safe_name

        existing_tree = fam_dir / "tree_rooted.nwk"
        if not existing_tree.exists():
            existing_tree = fam_dir / "tree.nwk"
        if args.skip_existing and existing_tree.exists() and (fam_dir / "annotations.tsv").exists():
            print(f"[{i+1}/{len(expanded_families)}] {family} — already exists, skipping")
            enzyme_name = family
            fam_parts   = family.split("|")
            for cid in central_seqs:
                cparts = cid.split("|")
                if len(cparts) >= 3 and cparts[0] == fam_parts[0] and cparts[1] == fam_parts[1]:
                    enzyme_name = re.sub(r"_\d+$", "", cparts[2])
                    break
            n_seqs = sum(1 for ln in open(fam_dir / "sequences.faa") if ln.startswith(">")) \
                if (fam_dir / "sequences.faa").exists() else 0
            families_built[family] = {"tree_file": str(existing_tree), "n_seqs": n_seqs,
                                      "enzyme_name": enzyme_name}
            continue

        print(f"[{i+1}/{len(expanded_families)}] {family} ({len(members)} proteins)")

        tree_seqs = collect_family_sequences(
            family, members, central_seqs, genome_seqs, mibig_seqs, mibig_hits,
        )
        print(f"  Collected {len(tree_seqs)} sequences for tree")

        if len(tree_seqs) < 4:
            log_skip(family, len(tree_seqs), "fewer than 4 sequences")
            families_skipped[family] = {"n_seqs": len(tree_seqs), "reason": "fewer than 4 sequences"}
            continue

        if args.max_seqs is not None and len(tree_seqs) > args.max_seqs:
            reason = f"too many sequences ({len(tree_seqs)} > {args.max_seqs})"
            log_skip(family, len(tree_seqs), reason)
            families_skipped[family] = {"n_seqs": len(tree_seqs), "reason": reason}
            continue

        tree_file = build_family_tree(family, tree_seqs, outdir, tools, threads=args.threads)

        if tree_file:
            # Identify EvoMining predictions (mirrors color_tree.pl)
            print("  Identifying EvoMining predictions via nw_clade...")
            evomining_preds = find_evomining_predictions(tree_file, tree_seqs, tools)
            for pred_id in evomining_preds:
                if pred_id in tree_seqs:
                    header, seq, _ = tree_seqs[pred_id]
                    tree_seqs[pred_id] = (header, seq, "evomining_prediction")
            print(f"  EvoMining predictions: {len(evomining_preds)} leaves marked green")

            enzyme_name = family
            fam_parts   = family.split("|")
            for cid in central_seqs:
                cparts = cid.split("|")
                if len(cparts) >= 3 and cparts[0] == fam_parts[0] and cparts[1] == fam_parts[1]:
                    enzyme_name = re.sub(r"_\d+$", "", cparts[2])
                    break

            # Written AFTER evomining_prediction reclassification
            render_tree_svg(family, tree_file, tree_seqs, outdir, tools)
            write_annotations(tree_seqs, outdir, family)
            if getattr(args, "itol", False):
                write_itol_colorstrip(tree_seqs, outdir, family)
            if getattr(args, "microreact", False):
                write_microreact_csv(tree_seqs, outdir, family)

            families_built[family] = {
                "tree_file":   tree_file,
                "n_seqs":      len(tree_seqs),
                "enzyme_name": enzyme_name,
            }
            print(f"  Tree built: {tree_file}")

    with open(skip_log_file, "w") as fh:
        fh.write("# EvoMining tree build — skipped families\n")
        fh.write(f"# Max sequences threshold: {args.max_seqs if args.max_seqs is not None else 'none (no limit)'}\n")
        fh.write(f"# Total skipped: {len(skip_log_lines)}\n#\n")
        fh.write("# status\tfamily\tsequences\treason\n")
        for line in skip_log_lines:
            fh.write(line + "\n")

    write_tree_index_html(families_built, families_skipped, outdir)

    print(f"\n{'='*60}")
    print(f"  Built {len(families_built)} trees")
    print(f"  Skipped {len(families_skipped)} families (see {skip_log_file})")
    print(f"  Output: {outdir}")
    print(f"{'='*60}")
    return outdir
