"""
evomining.cli
=============

The single `evomining` entry point.

There is no project file and no `init`. The current working directory is the
workspace: each command reads and writes at fixed conventional locations
beneath it (see evomining.workspace), so the pipeline chains without repeating
paths. Every default can still be overridden with the matching flag.

Typical session, run from your working folder:

    evomining generate-genome-db -i /vol/databases/.../genomes_annotated \\
                                 --list-dir ../0.clades/symb_clades
    evomining generate-enzyme-db --fasta core_symbiont_proteins_clean.faa
    evomining start   --mibig /path/to/MiBIG_DB.faa
    evomining analyze --mibig /path/to/MiBIG_DB.faa
    evomining trees   --mibig /path/to/MiBIG_DB.faa
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__, workspace
from .external import Tools, ToolError
from .workspace import WorkspaceError


def build_parser():
    ap = argparse.ArgumentParser(
        prog="evomining",
        description="EvoMining: enzyme family expansion analysis for BGC discovery.",
    )
    ap.add_argument("--version", action="version", version=f"evomining {__version__}")
    sub = ap.add_subparsers(dest="command", required=True, metavar="<command>")

    # -- check --
    p = sub.add_parser("check", help="verify every external tool is callable")
    p.add_argument("--tool-env", metavar="TOOL=PREFIX", action="append", default=[],
                   help="run TOOL from a separate conda env "
                        "(e.g. blastp=/vol/local/conda_envs/blast+)")

    # -- generate-genome-db --
    p = sub.add_parser("generate-genome-db", aliases=["genome-db"],
                       help="GenBank (.gbff) genomes -> EvoMining protein DB (no RAST)")
    p.add_argument("-i", "--input-dir", required=True,
                   help="directory of GenBank (.gbff/.gbk, optionally .gz) files, one per genome")
    p.add_argument("-l", "--lists", nargs="+", default=None,
                   help="restrict to these genome stems; .txt files, one stem per line")
    p.add_argument("--list-dir", default=None,
                   help="directory of .txt genome-stem lists (basename = clade)")
    p.add_argument("--names", default=None,
                   help="optional TSV overriding organism names (corason-py --names format)")
    p.add_argument("-o", "--output-dir", default=None,
                   help="output directory (default: ./evomining_db)")
    p.add_argument("--keep-pseudogenes", action="store_true",
                   help="keep pseudogene CDS (default: drop them; they are not functional "
                        "enzymes and can inflate copy counts). Use for continuity with older runs.")
    p.add_argument("-j", "--jobs", type=int, default=8,
                   help="number of parallel jobs to use (default: 8)")

    # -- generate-enzyme-db --
    p = sub.add_parser("generate-enzyme-db", aliases=["enzyme-db"],
                       help="build the Central (Enzyme) DB")
    g = p.add_mutually_exclusive_group(required=True)
    # General input modes (most users)
    g.add_argument("--fasta", metavar="FILE",
                   help="a Central DB FASTA already in canonical EvoMining format "
                        "(SUBSYSTEM|family_number|Function_querynumber|Organism); "
                        "headers are passed through verbatim")
    g.add_argument("--fasta-dir", metavar="DIR",
                   help="a directory of FASTA files, one file per enzyme family "
                        "(all sequences in a file are that family's seeds). Filename "
                        "encodes the family: SUBSYSTEM__enzyme.faa groups under a "
                        "subsystem, a bare enzyme.faa uses --custom-pathway-name")
    g.add_argument("--tsv", metavar="FILE",
                   help="a TSV (protein_id, family, enzyme_name, [subsystem]) whose "
                        "protein_ids are pulled from the genome FASTAs in --genomes-dir; "
                        "the optional 4th subsystem column groups families")
    # Specialized: cycad-transcriptome workflow
    g.add_argument("--transcripts", metavar="FILE",
                   help="SPECIALIZED (transcript workflow): a cleaned transcript FASTA "
                        "(TRINITY_ID|gene|host|description); each transcript becomes its own "
                        "single-seed custom|N family")
    # Options
    p.add_argument("--genomes-dir", default=None,
                   help="genome FASTA directory (required for --tsv mode)")
    p.add_argument("--custom-pathway-name", default="custom",
                   help="subsystem name used for families that have no subsystem of "
                        "their own (bare --fasta-dir filenames, --tsv rows without a "
                        "subsystem column, and all --transcripts families). Replaces the "
                        "'custom' in custom|N. Families that DO carry a subsystem keep it. "
                        "(default: custom)")
    p.add_argument("-o", "--output", default=None,
                   help="output Central DB FASTA (default: ./evomining_db/enzymes_db.faa)")

    # -- generate-antismash-db --
    p = sub.add_parser("generate-antismash-db", aliases=["antismash-db"],
                       help="build the OPTIONAL antiSMASH NP mapping (composite IDs)")
    p.add_argument("--antismash-dir", required=True,
                   help="antiSMASH results dir (per-genome dirs, optionally under clade subdirs)")
    p.add_argument("--genome-names", default=None,
                   help="genome_names.tsv from generate-genome-db "
                        "(default: ./evomining_db/genome_names.tsv)")
    p.add_argument("--genome-functions", default=None,
                   help="genome_functions.tsv from generate-genome-db "
                        "(default: alongside --genome-names)")
    p.add_argument("-o", "--output", default=None,
                   help="output mapping file (default: ./evomining_db/antismash_db.tsv)")

    # -- start --
    p = sub.add_parser("start", help="forward + reverse BLAST (Central DB <-> genome DB)")
    p.add_argument("-g", "--genomes", default=None,
                   help="genome protein FASTA (default: ./evomining_db/GENOMES.fasta)")
    p.add_argument("-c", "--central-db", default=None,
                   help="Central DB (default: ./evomining_db/enzymes_db.faa)")
    p.add_argument("-o", "--output-dir", default=None, help="default: ./evomining_out")
    p.add_argument("--threads", type=int, default=32, help="BLAST threads (default: 32)")
    p.add_argument("--evalue", type=float, default=1e-4, help="BLAST e-value (default: 1e-4)")
    p.add_argument("--chunk-size", type=int, default=100_000,
                   help="reverse BLAST query block size in proteins "
                        "(lower if blastp fails to create threads on a huge query; default: 100000)")

    # -- analyze --
    p = sub.add_parser("analyze", aliases=["analyse"],
                       help="copy counts, expansions, classification, MIBiG BLAST")
    p.add_argument("--blast-dir", default=None,
                   help="default: the single run under ./evomining_out")
    p.add_argument("--central-db", default=None)
    p.add_argument("--genome-names", default=None,
                   help="genome_names.tsv from generate-genome-db "
                        "(default: ./evomining_db/genome_names.tsv)")
    p.add_argument("--genome-functions", default=None,
                   help="genome_functions.tsv from generate-genome-db "
                        "(default: ./evomining_db/genome_functions.tsv)")
    p.add_argument("--genomes-fasta", default=None,
                   help="genome protein FASTA (default: ./evomining_db/GENOMES.fasta)")
    p.add_argument("--antismash", default=None)
    p.add_argument("--mibig", default=None, help="MIBiG protein FASTA for the NP BLAST")
    p.add_argument("-o", "--output-dir", default=None, help="default: ./evomining_analysis")
    p.add_argument("--threads", type=int, default=32)
    p.add_argument("--sd", type=float, default=1.0,
                   help="expansion threshold = mean + SD*stdev (default: 1.0)")
    p.add_argument("--ccm", default="multi", choices=["multi", "best"],
                help="Copy count method: multi (original, count protein in all matching families) or best (assign each protein to best-scoring family only)")
    # -- trees --
    p = sub.add_parser("trees", help="build colored phylogenetic trees per expanded family")
    p.add_argument("--classified", default=None,
                   help="default: ./evomining_analysis/classified_proteins.tsv")
    p.add_argument("--central-db", default=None)
    p.add_argument("--genomes-fasta", default=None,
                   help="genome protein FASTA (default: ./evomining_db/GENOMES.fasta)")
    p.add_argument("--mibig", default=None,
                   help="MIBiG FASTA; must match the one used by analyze")
    p.add_argument("--mibig-blast", default=None,
                   help="default: ./evomining_analysis/expanded_vs_mibig.blast")
    p.add_argument("--families", nargs="*", default=None,
                   help="only these families, or a .txt with one family per line")
    p.add_argument("--itol", action="store_true",
                   help="also write itol_colors.txt (iTOL TREE_COLORS) coloring each "
                        "leaf's branch/label by classification")
    p.add_argument("--microreact", action="store_true",
                   help="also write microreact.csv (id, classification, color) for "
                        "loading into microreact.org together with the .nwk tree")
    p.add_argument("--skip-existing", action="store_true")
    p.add_argument("--max-seqs", type=int, default=None,
                   help="skip families with more than this many sequences "
                        "(default: no limit)")
    p.add_argument("--threads", type=int, default=32)
    p.add_argument("-o", "--output-dir", default=None, help="default: ./evomining_trees")

    return ap


# ---------------------------------------------------------------------------
#  Commands
# ---------------------------------------------------------------------------

def cmd_check(args):
    tool_env = {}
    for item in args.tool_env:
        if "=" not in item:
            raise WorkspaceError(f"--tool-env expects TOOL=PREFIX, got: {item}")
        tool, prefix = item.split("=", 1)
        tool_env[tool] = prefix

    tools = Tools(tool_env)
    needed = ["makeblastdb", "blastp", "muscle", "mafft", "trimal", "FastTree",
              "nw_clade", "nw_labels", "nw_reroot", "nw_rename", "nw_display"]
    print("external tools:")
    missing = []
    for tool in needed:
        ok = tools.available(tool)
        where = tool_env.get(tool, "$PATH")
        print(f"  [{'ok' if ok else '--'}] {tool:<12} {where if ok else 'NOT FOUND'}")
        if not ok:
            missing.append(tool)
    if "FastTree" in missing and tools.available("fasttree"):
        missing.remove("FastTree")
    print()
    if missing:
        print(f"missing: {', '.join(missing)}")
        return 1
    print("all tools resolved")
    return 0


def _dir(override, name, flag, produced_by, what):
    return str(workspace.require_dir(
        workspace.resolve(override, name), flag=flag,
        produced_by=produced_by, what=what))


def _file(override, name, flag, produced_by, what):
    return str(workspace.require_file(
        workspace.resolve(override, name), flag=flag,
        produced_by=produced_by, what=what))


def _optional(override, name):
    """Override, else the conventional path only if it exists, else None."""
    if override:
        return override
    conv = workspace.default(name)
    return str(conv) if conv.exists() else None


def cmd_genome_db(args):
    from . import genomedb
    args.output_dir = args.output_dir or str(workspace.default("DB_DIR"))
    genomedb.run(args)
    return 0


def cmd_enzyme_db(args):
    from . import enzymedb
    args.output = args.output or str(workspace.default("ENZYME_DB"))
    # Ensure the parent directory exists (enzymes_db.faa lives under evomining_db/)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    if args.tsv and not args.genomes_dir:
        args.genomes_dir = _dir(args.genomes_dir, "GENOMES_DIR",
                                "--genomes-dir", "generate-genome-db", "GENOMES directory")
    enzymedb.run(args)
    return 0


def cmd_antismash_db(args):
    from . import antismashdb
    args.genome_names = _file(args.genome_names, "GENOME_NAMES",
                              "--genome-names", "generate-genome-db", "genome_names.tsv")
    # genome_functions.tsv sits beside genome_names.tsv in the DB dir; derive it from
    # the (possibly overridden) names path rather than assume a workspace key exists.
    gf = args.genome_functions or (Path(args.genome_names).parent / "genome_functions.tsv")
    args.genome_functions = str(workspace.require_file(
        Path(gf), flag="--genome-functions", produced_by="generate-genome-db",
        what="genome_functions.tsv"))
    args.output = args.output or str(workspace.default("ANTISMASH_DB"))
    antismashdb.run(args)
    return 0


def cmd_start(args):
    from . import start
    args.genomes = _file(args.genomes, "GENOMES_FASTA",
                         "--genomes", "generate-genome-db", "genome protein FASTA")
    args.central_db = _file(args.central_db, "ENZYME_DB",
                            "--central-db", "generate-enzyme-db", "Central (enzyme) DB")
    args.output_dir = args.output_dir or str(workspace.default("BLAST_ROOT"))
    start.run(args, tools=Tools())
    return 0


def cmd_analyze(args):
    from . import analysis
    args.blast_dir = str(workspace.find_blast_dir(args.blast_dir))
    args.central_db = _file(args.central_db, "ENZYME_DB",
                            "--central-db", "generate-enzyme-db", "Central (enzyme) DB")
    # Composite-ID artifacts from generate-genome-db (replaces --rast-ids / --genomes-dir).
    args.genome_names = _file(args.genome_names, "GENOME_NAMES",
                              "--genome-names", "generate-genome-db", "genome_names.tsv")
    # genome_functions.tsv sits beside genome_names.tsv in the DB dir; derive it from
    # the (possibly overridden) names path rather than assume a workspace key exists.
    gf = args.genome_functions or (Path(args.genome_names).parent / "genome_functions.tsv")
    args.genome_functions = str(workspace.require_file(
        Path(gf), flag="--genome-functions", produced_by="generate-genome-db",
        what="genome_functions.tsv"))
    args.genomes_fasta = _file(args.genomes_fasta, "GENOMES_FASTA",
                               "--genomes-fasta", "generate-genome-db",
                               "genome protein FASTA (GENOMES.fasta)")
    args.antismash = _optional(args.antismash, "ANTISMASH_DB")
    args.output_dir = args.output_dir or str(workspace.default("ANALYSIS_DIR"))
    analysis.run(args, tools=Tools())
    return 0


def cmd_trees(args):
    from . import trees
    classified = args.classified or (workspace.default("ANALYSIS_DIR") / "classified_proteins.tsv")
    args.classified = str(workspace.require_file(
        Path(classified), flag="--classified", produced_by="analyze",
        what="classified_proteins.tsv"))
    args.central_db = _file(args.central_db, "ENZYME_DB",
                            "--central-db", "generate-enzyme-db", "Central (enzyme) DB")
    args.genomes_fasta = _file(args.genomes_fasta, "GENOMES_FASTA",
                               "--genomes-fasta", "generate-genome-db",
                               "genome protein FASTA (GENOMES.fasta)")
    args.mibig_blast = _optional(args.mibig_blast, "MIBIG_BLAST")
    args.output_dir = args.output_dir or str(workspace.default("TREES_DIR"))
    trees.run(args, tools=Tools())
    return 0


DISPATCH = {
    "check": cmd_check,
    "generate-genome-db": cmd_genome_db, "genome-db": cmd_genome_db,
    "generate-enzyme-db": cmd_enzyme_db, "enzyme-db": cmd_enzyme_db,
    "generate-antismash-db": cmd_antismash_db, "antismash-db": cmd_antismash_db,
    "start": cmd_start,
    "analyze": cmd_analyze, "analyse": cmd_analyze,
    "trees": cmd_trees,
}


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return DISPATCH[args.command](args)
    except (WorkspaceError, ToolError) as exc:
        print(f"\nerror: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
