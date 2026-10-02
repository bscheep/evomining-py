#!/usr/bin/env python3
"""
evomining.enzymedb

Build a custom EvoMining Central (Enzyme) DB. The Central DB defines which
enzyme families EvoMining tracks across your genome database; for each family
you provide 2-5 seed sequences.

Central DB header format (canonical EvoMining):

    SUBSYSTEM|family_number|Function_querynumber|Organism_or_comment

  SUBSYSTEM          metabolic subsystem; a subsystem may hold many families.
  family_number      defines a column on the heatplot. A family may have many
                     query seeds; SUBSYSTEM|family_number is the grouping key.
  Function_query...   the family's function, plus a unique consecutive query id.
  Organism_or_comment organism name or free comment.

  Example:  >3PGA_AMINOACIDS|1|Phosphoglyceratedehydrogenase_1|Cglu

Two seeds of the SAME family share SUBSYSTEM|family_number and differ only in
the querynumber suffix / organism, e.g.

    >3PGA_AMINOACIDS|1|Phosphoglyceratedehydrogenase_1|Cglu
    >3PGA_AMINOACIDS|1|Phosphoglyceratedehydrogenase_2|Ecoli

This is exactly the key analysis.parse_central_header reads, so a DB built here
is consumed downstream without translation.

Input modes:
  --fasta        A FASTA already in canonical format (headers passed through
                 verbatim; your family numbering and organism strings are kept).
  --fasta-dir    A directory of FASTAs, one file per family. All sequences in a
                 file are that family's seeds. Filename encodes the family, with
                 an optional subsystem prefix: `SUBSYSTEM__enzyme.faa` groups the
                 family under SUBSYSTEM; a bare `enzyme.faa` uses the default
                 subsystem (--custom-pathway-name).
  --tsv          A TSV (protein_id, family, enzyme_name, [subsystem]) + a genome
                 FASTA dir. The optional 4th subsystem column groups families;
                 rows sharing a family accumulate as multiple seeds.
  --transcripts  SPECIALIZED (cycad-transcriptome workflow): a cleaned transcript
                 FASTA whose headers are TRINITY_ID|gene|host|description. Each
                 transcript becomes its own single-seed custom|N family.

Usage:
  evomining generate-enzyme-db --fasta central_db.faa
  evomining generate-enzyme-db --fasta-dir enzymes/ -o my_central_db.faa
  evomining generate-enzyme-db --transcripts core_symbiont_proteins_clean.faa
"""

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path


def read_fasta(path):
    """Read FASTA -> list of (header, sequence)."""
    entries = []
    header = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    entries.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header:
        entries.append((header, "".join(seq)))
    return entries


def sanitize_name(name):
    """Sanitize a name for use in FASTA headers."""
    name = re.sub(r"[^a-zA-Z0-9_\-.]", "_", name)
    name = re.sub(r"-_", "-", name)   # remove underscore after hyphen
    name = re.sub(r"_-", "-", name)   # remove underscore before hyphen
    while "__" in name:
        name = name.replace("__", "_")
    return name.strip("_")


# ---------------------------------------------------------------------------
#  Mode: canonical FASTA (the default for --fasta)
# ---------------------------------------------------------------------------

def from_canonical_fasta(fasta_file):
    """
    Read a Central DB FASTA that is ALREADY in canonical EvoMining format:

        >SUBSYSTEM|family_number|Function_querynumber|Organism_or_comment

    Field 0 (subsystem) + field 1 (family number) together define the enzyme
    family, matching analysis.parse_central_header. Field 2 is the function
    label, field 3 the organism/comment (optional). Headers are preserved
    verbatim -- the user's own family numbering and organism strings are
    respected, never reassigned.

    Returns (entries, summary):
      entries -- list of (header, seq) with the ORIGINAL header preserved
      summary -- dict of subsystem / family / seed counts for reporting

    Hard-fails on any header with fewer than 3 pipe fields rather than
    force-fitting it, so a non-canonical file is caught early with a clear
    message (use --transcripts for cleaned transcript headers).
    """
    entries = read_fasta(fasta_file)
    if not entries:
        raise SystemExit(f"ERROR: no sequences in {fasta_file}")

    family_seed_counts = defaultdict(int)
    subsystems = set()
    bad = []
    whitespace = []
    for header, seq in entries:
        fields = header.split("|")
        if len(fields) < 3:
            bad.append(header)
            continue
        subsystems.add(fields[0])
        family_seed_counts[f"{fields[0]}|{fields[1]}"] += 1
        # BLAST truncates sequence IDs at the first whitespace, so a space in a
        # canonical header would silently drop the family/query fields after it.
        if any(ch.isspace() for ch in header):
            whitespace.append(header)

    if bad:
        preview = "\n    ".join(bad[:5])
        more = f"\n    ... and {len(bad) - 5} more" if len(bad) > 5 else ""
        raise SystemExit(
            f"ERROR: {len(bad)} header(s) are not in canonical EvoMining format.\n"
            f"  Expected: SUBSYSTEM|family_number|Function_querynumber|Organism\n"
            f"  Example:  3PGA_AMINOACIDS|1|Phosphoglyceratedehydrogenase_1|Cglu\n"
            f"  Offending:\n    {preview}{more}\n"
            f"  (For cleaned transcript headers TRINITY_ID|gene|host|description, "
            f"use --transcripts instead.)"
        )

    if whitespace:
        print(f"  WARNING: {len(whitespace)} header(s) contain whitespace; BLAST will "
              f"truncate the sequence ID at the first space, dropping later fields. "
              f"Replace spaces with underscores.", file=sys.stderr)

    summary = {
        "n_subsystems": len(subsystems),
        "n_families": len(family_seed_counts),
        "n_seeds": len(entries),
    }
    return entries, summary


def write_canonical_db(entries, output_file):
    """Write canonical entries verbatim (headers untouched, sequence wrapped at 60)."""
    output = Path(output_file)
    with open(output, "w") as fh:
        for header, seq in entries:
            fh.write(f">{header}\n")
            for j in range(0, len(seq), 60):
                fh.write(seq[j:j + 60] + "\n")
    return output


# ---------------------------------------------------------------------------
#  Mode: directory of per-family FASTAs
# ---------------------------------------------------------------------------

def from_fasta_dir(fasta_dir, pathway_name="custom"):
    """
    Build Central DB from a directory of FASTA files, one file per enzyme family.
    All sequences in a file become that family's seeds.

    The filename encodes the family and its function, optionally with a
    subsystem prefix:

      SUBSYSTEM__enzyme.faa   ->  family key "SUBSYSTEM|<n>"  (subsystem taken
                                  from the part before the '__' separator;
                                  `enzyme` becomes the Function field shared
                                  by every seed in the family)
      enzyme.faa              ->  family key "enzyme"         (no subsystem; falls
                                  back to --custom-pathway-name at write time)

    e.g. `3PGA_AMINOACIDS__cysteine_synthase.faa` groups under subsystem
    3PGA_AMINOACIDS with Function "cysteine_synthase", while
    `cysteine_synthase.faa` gets the default pathway name.

    Each sequence's own header is carried through UNCHANGED and becomes that
    seed's Organism/comment field -- it is never parsed, split, or discarded.
    This is what lets a composite `<genome_stem>__<locus_tag>` header end up
    as the seed's organism identity downstream (trees.py keys its Central DB
    seed lookup on the full header, and uses this field to label tree leaves).
    """
    families = {}
    fasta_dir = Path(fasta_dir)

    for fasta_file in sorted(fasta_dir.glob("*.faa")) + sorted(fasta_dir.glob("*.fasta")) + sorted(fasta_dir.glob("*.fa")):
        stem = fasta_file.stem
        entries = read_fasta(str(fasta_file))
        if not entries:
            continue
        if "__" in stem:
            subsystem, enzyme = stem.split("__", 1)
            subsystem = sanitize_name(subsystem)
            enzyme = sanitize_name(enzyme)
            # Key on "subsystem|enzyme" so write_central_db keeps the subsystem
            # prefix; the enzyme part becomes the shared Function label.
            family_key = f"{subsystem}|{enzyme}"
        else:
            enzyme = sanitize_name(stem)
            family_key = enzyme

        seen_headers = set()
        for header, _ in entries:
            if header in seen_headers:
                print(f"  WARNING: duplicate sequence header {header!r} in "
                      f"{fasta_file.name} -- these seeds will collide "
                      f"downstream (no per-seed query number is added in "
                      f"this mode; uniqueness relies on distinct headers).",
                      file=sys.stderr)
            seen_headers.add(header)

        # (function_label, organism_field, seq) -- organism_field is the
        # seed's own header, untouched, relied on for per-seed uniqueness.
        families[family_key] = [(enzyme, header, seq) for header, seq in entries]

    return families


# ---------------------------------------------------------------------------
#  Mode: cleaned transcript FASTA -> canonical custom|N families
# ---------------------------------------------------------------------------

def from_transcript_fasta(fasta_file):
    """
    Build Central DB seeds from a cleaned transcript FASTA.

    Handles cleaned transcript headers like:
      >TRINITY_DN101449_c0_g1_i1|SYT|Cycas|Threonine--tRNA_ligase

    Each transcript becomes its own single-seed family (keyed by the Trinity
    id); write_central_db then numbers them custom|1, custom|2, ... The seed
    label is built as gene_description then suffixed with the host, and a
    literal "unknown" description is dropped. Falls back to family|seed or
    whitespace-delimited headers for non-transcript inputs.
    """
    families = defaultdict(list)
    entries = read_fasta(fasta_file)

    for header, seq in entries:
        if "|" in header:
            parts = header.split("|")
            if len(parts) >= 4:
                # Cleaned transcript format: TRINITY_ID|gene_name|host|description
                trinity_id = parts[0]
                gene_name = parts[1]
                host = parts[2]
                description = parts[3] if len(parts) > 3 else ""
                if description and description != "unknown":
                    enzyme_name = f"{gene_name}_{description}"
                else:
                    enzyme_name = gene_name
                seed_name = f"{sanitize_name(enzyme_name)}_{host}"
                families[trinity_id].append((seed_name, seq))
            elif len(parts) >= 2:
                # Simple family|seed format
                family = sanitize_name(parts[0])
                seed_name = parts[1] if len(parts) > 1 else header
                families[family].append((seed_name, seq))
        else:
            parts = header.split(None, 1)
            family = sanitize_name(parts[0])
            seed_name = parts[1] if len(parts) > 1 else header
            families[family].append((seed_name, seq))

    return {fam: [(name, seq) for name, seq in seeds] for fam, seeds in families.items()}


# ---------------------------------------------------------------------------
#  Mode: TSV mapping + genome FASTA files
# ---------------------------------------------------------------------------

def from_tsv_mapping(tsv_file, genomes_dir):
    """
    Build Central DB from a TSV file mapping proteins to families.

    TSV columns (tab-separated, with a header row):
        protein_id  family  enzyme_name  [subsystem]

    The 4th `subsystem` column is optional. When present, the family key becomes
    "subsystem|family" so write_central_db keeps the subsystem prefix and numbers
    families within it; when absent, the family is plain "family" and falls back
    to --custom-pathway-name. Multiple rows sharing the same family (and subsystem)
    accumulate as multiple seeds of that family.

    protein_id should match a header in the genome .faa files.
    """
    families = defaultdict(list)
    genomes_dir = Path(genomes_dir)

    # Load all genome sequences into memory (indexed by leading id token)
    all_seqs = {}
    for faa in sorted(genomes_dir.glob("*.faa")):
        for header, seq in read_fasta(str(faa)):
            seq_id = header.split()[0]
            all_seqs[seq_id] = (header, seq)

    with open(tsv_file) as fh:
        header_line = fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            protein_id = parts[0]
            family = sanitize_name(parts[1])
            enzyme_name = parts[2]
            subsystem = sanitize_name(parts[3]) if len(parts) >= 4 and parts[3].strip() else None

            # Key on "subsystem|family" if a subsystem was given, else plain family.
            family_key = f"{subsystem}|{family}" if subsystem else family

            seq_entry = all_seqs.get(protein_id) or all_seqs.get(f"fig|{protein_id}")
            if seq_entry:
                _, seq = seq_entry
                families[family_key].append((enzyme_name, seq))
            else:
                print(f"  WARNING: protein not found: {protein_id}", file=sys.stderr)

    return dict(families)


# ---------------------------------------------------------------------------
#  Write Central DB (for the family-dict modes: transcripts / dir / tsv / extract)
# ---------------------------------------------------------------------------

def write_central_db(families, output_file, pathway_name="custom"):
    """
    Write EvoMining-compatible Central DB FASTA.

    Format: >SUBSYSTEM|family_number|Function|Organism

    Seeds from --fasta-dir carry a real per-seed Organism (the sequence's own
    header, untouched) and share one Function label per family (from the
    filename). Seeds from --tsv / --transcripts keep the legacy behaviour:
    Function gets a unique trailing query number and Organism falls back to
    the literal pathway_name (neither mode tracks a real organism per seed).

    The family_number must be an integer (it defines a heatplot column, and
    analysis.parse_central_header keys the family on SUBSYSTEM|family_number).

    Family key handling:
      "SUBSYSTEM|enzyme"  -> subsystem = SUBSYSTEM, family numbered per subsystem
      "enzyme"            -> subsystem = pathway_name (--custom-pathway-name),
                             family numbered globally

    So `3PGA_AMINOACIDS|cysteine_synthase` becomes `3PGA_AMINOACIDS|1|...`,
    a second family in the same subsystem becomes `3PGA_AMINOACIDS|2|...`, and a
    plain `citrate_synthase` becomes `<pathway_name>|<n>|...`.

    (Canonical --fasta input does NOT go through here; it is written verbatim by
    write_canonical_db so the user's own numbering/organism fields are preserved.)
    """
    output = Path(output_file)
    total_seqs = 0

    # Assign a per-subsystem family number counter.
    subsystem_counter = defaultdict(int)
    global_counter = 0

    with open(output, "w") as fh:
        for family_key in sorted(families.keys()):
            seeds = families[family_key]

            if "|" in family_key:
                subsystem, enzyme_label = family_key.split("|", 1)
            else:
                subsystem = pathway_name
                enzyme_label = family_key

            subsystem_counter[subsystem] += 1
            family_num = subsystem_counter[subsystem]

            for i, seed in enumerate(seeds):
                if isinstance(seed, tuple) and len(seed) == 3:
                    # fasta-dir mode: function is shared by the whole family
                    # (from the filename); organism is the seed's own header,
                    # already unique, so no query-number suffix is needed.
                    seed_name, organism, seq = seed
                    func_label = sanitize_name(seed_name)
                    organism_clean = sanitize_name(organism)
                elif isinstance(seed, tuple):
                    # tsv / transcripts mode: unchanged legacy behaviour.
                    seed_name, seq = seed
                    func_label = f"{sanitize_name(seed_name)}_{i+1}"
                    organism_clean = pathway_name
                else:
                    seq = seed
                    func_label = f"seed_{i+1}"
                    organism_clean = pathway_name

                header = f">{subsystem}|{family_num}|{func_label}|{organism_clean}"
                fh.write(f"{header}\n")
                for j in range(0, len(seq), 60):
                    fh.write(seq[j:j+60] + "\n")
                total_seqs += 1

    n_families = sum(subsystem_counter.values())
    return n_families, total_seqs


def run(args):
    """Build the Central (Enzyme) DB from canonical seeds, transcripts, or references."""
    print("Building EvoMining Enzyme DB...")
    n_subsystems = None

    if args.fasta:
        # Default: the file is already in canonical EvoMining format.
        print(f"  Mode: canonical EvoMining FASTA ({args.fasta})")
        entries, summary = from_canonical_fasta(args.fasta)
        write_canonical_db(entries, args.output)
        n_families = summary["n_families"]
        n_seqs = summary["n_seeds"]
        n_subsystems = summary["n_subsystems"]

    elif getattr(args, "transcripts", None):
        print(f"  Mode: cleaned transcript headers -> canonical ({args.transcripts})")
        raw = from_transcript_fasta(args.transcripts)
        families = {}
        for fam, seeds in raw.items():
            if seeds and isinstance(seeds[0], tuple):
                families[fam] = seeds
            else:
                families[fam] = [(f"seed_{i+1}", sq) for i, sq in enumerate(seeds)]
        if not families:
            raise SystemExit("ERROR: no enzyme families found")
        n_families, n_seqs = write_central_db(families, args.output, args.custom_pathway_name)

    elif args.fasta_dir:
        print(f"  Mode: directory of per-family FASTAs ({args.fasta_dir})")
        families = from_fasta_dir(args.fasta_dir, args.custom_pathway_name)
        if not families:
            raise SystemExit("ERROR: no enzyme families found")
        n_families, n_seqs = write_central_db(families, args.output, args.custom_pathway_name)

    elif args.tsv:
        if not args.genomes_dir:
            raise SystemExit("ERROR: --genomes-dir required for --tsv mode")
        print(f"  Mode: TSV mapping ({args.tsv})")
        families = from_tsv_mapping(args.tsv, args.genomes_dir)
        if not families:
            raise SystemExit("ERROR: no enzyme families found")
        n_families, n_seqs = write_central_db(families, args.output, args.custom_pathway_name)

    else:
        raise SystemExit("ERROR: pick one input mode "
                         "(--fasta / --transcripts / --fasta-dir / --tsv)")

    print(f"\n{'='*60}")
    print(f"  Enzyme DB built")
    if n_subsystems is not None:
        print(f"  Subsystems:      {n_subsystems}")
    print(f"  Enzyme families: {n_families}")
    print(f"  Total seeds:     {n_seqs}")
    print(f"  Output:          {args.output}")
    print(f"{'='*60}")
    return Path(args.output)
