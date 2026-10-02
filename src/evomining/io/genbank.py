"""GenBank flat-file loader.

Adapted from corason-py (miguel-mx/corason-py), src/corason/io/genbank.py. The
parsing logic is reused essentially unchanged so EvoMining and CORASON identify
and translate genes identically; the only change is that ``_finalize_contigs``
now lives in EvoMining's own ``io.models``.

Defensive about what differs between a RAST dump, an NCBI RefSeq download, a
Prokka run, and a fungal annotation:

* **Identifiers** fall back ``locus_tag`` -> ``protein_id`` -> ``gene`` ->
  synthesized. No integer is ever parsed out of an identifier.
* **Isoforms** -- several CDS sharing one ``locus_tag`` -- collapse to one gene.
* **Compound locations** (``join(...)``, spliced CDS) keep their exon extents.
* **Translations** come from ``/translation`` when present, else are derived from
  the contig sequence honoring ``/transl_table`` and ``/codon_start``.
* **Pseudogenes** are skipped rather than translated into nonsense.
* **Compressed input** (``.gbff.gz``, ``.gbk.gz``, ...) is read transparently.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from .models import Gene, Genome, finalize_contigs

logger = logging.getLogger(__name__)

GENBANK_SUFFIXES = frozenset({".gbk", ".gb", ".gbff", ".genbank"})

_DEFAULT_TRANSL_TABLE = 11  # bacterial / archaeal / plant plastid


def _plain(path: Path) -> Path:
    """``path`` with a trailing ``.gz`` removed: ``Foo.gbff.gz`` -> ``Foo.gbff``."""
    return path.with_suffix("") if path.suffix.lower() == ".gz" else path


def is_genbank(path: Path) -> bool:
    return _plain(path).suffix.lower() in GENBANK_SUFFIXES


def genbank_suffix(path: Path) -> str:
    """The GenBank suffix, ignoring a trailing ``.gz`` (``.gbff`` for ``x.gbff.gz``)."""
    return _plain(path).suffix.lower()


def genome_stem(path: Path) -> str:
    """The genome stem of a (possibly gzipped) GenBank file.

    ``Path.stem`` alone is wrong for compressed input -- ``Foo.gbff.gz`` would give
    ``Foo.gbff``. The stem namespaces every protein ID and has to match ``--lists``
    entries and antiSMASH directory names, so getting it wrong breaks those joins
    silently rather than loudly. Always go through this helper.
    """
    return _plain(path).stem


def _parse_records(path: Path):
    """Yield the GenBank records of ``path``, decompressing ``.gz`` transparently."""
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt") as handle:
        yield from SeqIO.parse(handle, "genbank")


def load_genbank(path: Path, name: str | None = None, keep_pseudo: bool = False,
                 stem: str | None = None) -> Genome:
    """Load every record of a GenBank file into one :class:`Genome`.

    ``stem`` overrides the genome id (which defaults to the filename stem). It is
    set when genomes are laid out one-per-folder, so the id becomes the folder
    name rather than a possibly-generic inner filename.
    """
    gid = stem or genome_stem(path)
    genome_name = name
    accessions: list[str] = []

    by_contig: dict[str, list[Gene]] = {}
    contig_lengths: dict[str, int | None] = {}
    counts = _SkipCounts(keep_pseudo)

    for record in _parse_records(path):
        contig_id = record.id or record.name
        if not contig_id or contig_id == "<unknown id>":
            contig_id = f"{gid}_{len(by_contig) + 1}"

        if record.id and record.id not in accessions:
            accessions.append(record.id)
        if not genome_name:
            genome_name = record.annotations.get("organism", "") or ""

        contig_lengths[contig_id] = len(record.seq) if _has_sequence(record.seq) else None
        genes: list[Gene] = []

        for feature in record.features:
            if feature.type != "CDS":
                continue
            gene = _build_gene(feature, record, contig_id, len(genes), counts, keep_pseudo)
            if gene is not None:
                genes.append(gene)

        genes = _collapse_isoforms(genes, counts)
        if genes:
            by_contig[contig_id] = genes

    genome_name = genome_name or genome_stem(path)

    counts.report(path)

    genome = Genome(id=gid, name=genome_name, accessions=accessions, source_path=path)
    
    finalize_contigs(genome, by_contig)

    for contig_id, length in contig_lengths.items():
        if length is not None and contig_id in genome.contigs:
            genome.contigs[contig_id].length = length

    if genome.n_genes == 0:
        raise ValueError(f"{path}: no usable CDS features found")

    return genome


class _SkipCounts:
    """Tally of features dropped during parsing, reported once rather than per-feature."""

    def __init__(self, keep_pseudo: bool = False) -> None:
        self.keep_pseudo = keep_pseudo
        self.pseudo = 0
        self.untranslatable = 0
        self.ambiguous = 0
        self.no_strand = 0
        self.isoforms = 0

    def report(self, path: Path) -> None:
        if self.pseudo:
            verb = "kept" if self.keep_pseudo else "skipped"
            logger.warning("%s: %s %d pseudogene CDS", path.name, verb, self.pseudo)
        if self.isoforms:
            logger.warning(
                "%s: collapsed %d CDS sharing a locus_tag; kept the longest product",
                path.name, self.isoforms,
            )
        if self.untranslatable:
            logger.warning(
                "%s: skipped %d CDS with no /translation and no derivable sequence",
                path.name, self.untranslatable,
            )
        if self.ambiguous:
            logger.warning(
                "%s: skipped %d CDS lying entirely in an assembly gap (all-N)",
                path.name, self.ambiguous,
            )
        if self.no_strand:
            logger.warning(
                "%s: %d CDS had no strand; assumed forward", path.name, self.no_strand
            )


def _collapse_isoforms(genes: list[Gene], counts: _SkipCounts) -> list[Gene]:
    """One gene per identifier. Several CDS may legitimately share a ``locus_tag``
    (alternative start codons, programmed ribosomal frameshifts). The longest
    product wins, with coordinates and identifier breaking ties, so the choice
    never depends on feature order in the file."""
    best: dict[str, Gene] = {}
    order: list[str] = []

    for gene in genes:
        held = best.get(gene.id)
        if held is None:
            best[gene.id] = gene
            order.append(gene.id)
            continue
        counts.isoforms += 1
        if _isoform_key(gene) < _isoform_key(held):
            best[gene.id] = gene

    return [best[gene_id] for gene_id in order]


def _isoform_key(gene: Gene) -> tuple[int, int, int, str]:
    return (-len(gene.translation), gene.start, gene.end, gene.translation)


def _has_sequence(seq: Seq | None) -> bool:
    """Whether a record carries real residues, as opposed to a CONTIG stub."""
    if seq is None or len(seq) == 0:
        return False
    try:
        bytes(seq)
    except Exception:  # Bio.Seq raises UndefinedSequenceError for CONTIG-only records
        return False
    return True


def _build_gene(feature, record, contig_id: str, ordinal: int, counts: _SkipCounts,
                keep_pseudo: bool = False):
    is_pseudo = "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers
    if is_pseudo:
        counts.pseudo += 1
        if not keep_pseudo:
            return None
        # kept for continuity: a pseudogene often has no /translation, so it may
        # still be dropped below by _translation() if nothing is derivable.

    translation = _translation(feature, record, counts)
    if not translation:
        return None

    location = feature.location
    parts = [(int(p.start), int(p.end)) for p in location.parts]
    is_compound = len(parts) > 1

    strand = location.strand
    if strand is None:
        counts.no_strand += 1
        strand = 1

    return Gene(
        id=_gene_id(feature, contig_id, ordinal),
        contig=contig_id,
        index=-1,  # assigned by finalize_contigs, after sorting
        start=int(location.start),
        end=int(location.end),
        strand=1 if strand >= 0 else -1,
        product=_first(feature, "product") or "hypothetical protein",
        translation=translation,
        is_compound=is_compound,
        parts=parts if is_compound else [],
    )


def _gene_id(feature, contig_id: str, ordinal: int) -> str:
    """Stable identifier, preferring the most specific qualifier available."""
    for key in ("locus_tag", "protein_id", "gene"):
        value = _first(feature, key)
        if value:
            return value
    return f"{contig_id}_cds_{ordinal + 1}"


def _first(feature, key: str) -> str | None:
    values = feature.qualifiers.get(key)
    return values[0].strip() if values else None


def _translation(feature, record, counts: _SkipCounts) -> str | None:
    """The protein sequence, from ``/translation`` or derived from the contig,
    honoring ``/transl_table`` and ``/codon_start``."""
    stated = _first(feature, "translation")
    if stated:
        return stated.replace(" ", "").rstrip("*")

    if not _has_sequence(record.seq):
        counts.untranslatable += 1
        return None

    table = _first(feature, "transl_table") or str(_DEFAULT_TRANSL_TABLE)
    try:
        table_id = int(table)
        nucleotides = feature.extract(record.seq)
        offset = int(_first(feature, "codon_start") or 1) - 1
        nucleotides = nucleotides[offset:]
    except Exception:
        counts.untranslatable += 1
        return None

    try:
        protein = nucleotides.translate(table=table_id, cds=True)
    except Exception:
        whole = nucleotides[: len(nucleotides) // 3 * 3]
        try:
            protein = whole.translate(table=table_id, to_stop=True)
        except Exception:
            counts.untranslatable += 1
            return None

    protein = str(protein).rstrip("*")
    if not protein:
        counts.untranslatable += 1
        return None

    if set(protein) <= {"X"}:
        counts.ambiguous += 1
        return None

    return protein
