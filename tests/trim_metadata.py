from pathlib import Path

BLAST_DIR = Path("data/blast")
METADATA_DIR = Path("data/metadata")

referenced = set()
for name, col in [("thereand.blast", 1), ("backagain.blast", 0)]:
    with open(BLAST_DIR / name) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 12:
                value = fields[col]
                if "__" in value:
                    referenced.add(value)

print(f"{len(referenced)} referenced proteins found")

for name in ("genome_names.tsv", "genome_functions.tsv"):
    path = METADATA_DIR / name
    with open(path) as fh:
        lines = fh.readlines()
    header, rows = lines[0], lines[1:]
    kept = [header] + [r for r in rows if r.split("\t", 1)[0] in referenced]
    with open(path, "w") as fh:
        fh.writelines(kept)
    print(f"{name}: {len(rows)} -> {len(kept) - 1} rows")
