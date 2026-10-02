#!/usr/bin/env python3
"""
evomining.genomedb
==================

Build EvoMining's genome database from GenBank (.gbff/.gbk) input, given either as
a flat folder of GenBank files (stem = filename) or as a folder of per-genome
sub-folders each holding that genome's annotation (stem = folder name, e.g. Bakta
output). Both layouts are accepted transparently; see resolve_genome_inputs.

This replaces the former RAST-table builder. There is no more 13-column .txt, no
Corason_Rast.IDs, and no 666666.<id>.peg.<N> scheme. Genomes are parsed straight
from GenBank by the vendored corason-py parser (evomining.io), and EvoMining's
downstream steps consume three flat artifacts written here:

    GENOMES.fasta          protein FASTA, one record per CDS, header = composite ID
    genome_functions.tsv   composite_id <TAB> product
    genome_names.tsv       genome_id <TAB> organism_display_name

The composite protein ID is:

    <genome_stem>__<locus_tag>

i.e. the GenBank filename stem, a double underscore, and the gene's stable
identifier (locus_tag / protein_id / gene / synthesized, per the parser). The
stem namespaces the locus so that genomes which reuse locus tags (Bakta/Prokka
restart numbering per genome) do not collide when flattened into one FASTA. This
matches corason-py's own <genome>__<locus_tag> convention for its .gbk/cluster
names. The stem is used (not the organism name) so IDs stay stable regardless of
what else is in the batch, since organism names may be disambiguated at load time.

Optional --names TSV overrides organism display names for genomes with poor
GenBank metadata; it uses the same format as corason-py's --names/--rast-ids.

Clade / iTOL metadata is intentionally NOT produced here: clade assignment is a
study-specific overlay, not part of EvoMining. Build any clade colorstrips with a
separate tool from the genome lists.
"""

from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from collections import deque
from itertools import islice
from pathlib import Path

from .utils import setup_logging
from .io.genbank import load_genbank
from .io.loader import disambiguate_names, resolve_genome_inputs
from .io.names import load_names



ID_SEP = "__"


def composite_id(genome_stem: str, gene_id: str) -> str:
    return f"{genome_stem}{ID_SEP}{gene_id}"

def _patch_genome_names(names_path: Path, renamed: dict[str, str]):
    """Rewrite the genome_names.tsv for genomes that have been renamed."""
    tmp_path = names_path.with_suffix(names_path.suffix + ".tmp")
    with open(names_path) as src, open(tmp_path, "w") as dst:
        for line in src:
            stem = line.split("\t", 1)[0]
            if stem in renamed:
                dst.write(f"{stem}\t{renamed[stem]}\n")
            else:
                dst.write(line)
    tmp_path.replace(names_path)


def _wanted_stems(args):
    """Set of genome stems to include from --lists / --list-dir, or None for all."""
    stems = set()
    have_lists = False

    if args.lists:
        have_lists = True
        for lst in args.lists:
            stems |= _read_list(Path(lst))
    if args.list_dir:
        have_lists = True
        for lst in sorted(Path(args.list_dir).glob("*.txt")):
            stems |= _read_list(lst)

    return stems if have_lists else None


def _read_list(path):
    out = set()
    with open(path) as fh:
        for line in fh:
            name = line.strip()
            if name:
                out.add(name)
    return out


def _parse_all(pairs, jobs, name_overrides, keep_pseudo):
    """Yield one Genome per (path, stem) pair, in 'pairs' order."""
    if jobs == 1:
        for path, stem in pairs:
            yield load_genbank(path, name=name_overrides.get(stem), keep_pseudo=keep_pseudo, stem=stem)
        return
    
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        def submit(path, stem):
            return pool.submit(load_genbank, path, name=name_overrides.get(stem), keep_pseudo=keep_pseudo, stem=stem)
        
        it = iter(pairs)
        pending = deque(submit(p, s) for p, s in islice(it, jobs * 2))
        while pending:
            genome = pending.popleft().result()
            nxt = next(it, None)
            if nxt is not None:
                pending.append(submit(*nxt))
            yield genome


def run(args):
    """Parse GenBank genomes and write EvoMining's flat genome artifacts."""
    input_dir = Path(args.input_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(outdir / "evomining.log", name="evomining")
    try:
        pairs = resolve_genome_inputs(input_dir)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}")

    wanted = _wanted_stems(args)
    if wanted is not None:
        pairs = [(p, s) for (p, s) in pairs if s in wanted]
        missing = wanted - {s for _, s in pairs}
        for stem in sorted(missing):
            print(f"  WARNING: no GenBank file for: {stem}", file=sys.stderr)

    if not pairs:
        raise SystemExit("ERROR: no genomes selected")

    keep_pseudo = getattr(args, 'keep_pseudogenes', False)
    name_overrides = load_names(Path(args.names)) if args.names else {}

    fasta_path = outdir / "GENOMES.fasta"
    func_path = outdir / "genome_functions.tsv"
    names_path = outdir / "genome_names.tsv"

    logger.info(f"Loading {len(pairs)} GenBank genome(s)...")
    metas = []
    total_genes = 0
    with open(fasta_path, "w") as fa, \
         open(func_path, "w") as fn, \
         open(names_path, "w") as nm:
        fn.write("protein_id\tfunction\n")
        nm.write("genome_id\tgenome_name\n")
        for genome in _parse_all(pairs, getattr(args, 'jobs', 8), name_overrides, keep_pseudo):
            nm.write(f"{genome.metadata.id}\t{genome.metadata.name}\n")
            gene_count = 0
            for gene in genome.genes():
                cid = composite_id(genome.metadata.id, gene.id)
                fa.write(f">{cid}\n")
                for i in range(0, len(gene.translation), 60):
                    fa.write(gene.translation[i:i + 60] + "\n")
                fn.write(f"{cid}\t{gene.product}\n")
                gene_count += 1
            total_genes += gene_count
            logger.info(f"  [{genome.metadata.id}]  {genome.metadata.name:40s}  {gene_count:>5} proteins  "
                        f"({len(genome.contigs)} contigs)")
            metas.append(genome.metadata)

    before = {m.id: m.name for m in metas}
    disambiguate_names(metas)
    renamed = {m.id: m.name for m in metas if m.name != before[m.id]}
    if renamed:
        logger.info(f"Fixing up {len(renamed)} duplicate organism name(s) in {names_path.name}...")
        _patch_genome_names(names_path, renamed)


    logger.info(f"""
{'=' * 60}
  Done!
  Genomes:           {len(metas)}
  Total proteins:    {total_genes}

  Protein FASTA:     {fasta_path}
  Function map:      {func_path}
  Name map:          {names_path}
{'=' * 60}
""")
    return outdir
