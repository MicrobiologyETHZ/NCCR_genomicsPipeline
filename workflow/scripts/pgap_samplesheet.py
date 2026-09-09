"""
Build a PGAP samplesheet from a directory of genome FASTAs.

PGAP needs a per-genome taxon (`-s`), which can't be inferred from a filename,
so this writes an editable CSV rather than driving PGAP straight off a
directory scan: run this once, then hand-edit the `taxon` column for any
genome that needs a different value than the `--taxon` default.

Discovery reuses phage_inputs.resolve_assemblies for its collision-safe naming
and FASTA-pattern matching, with `recursive: True` so subdirectories are
included. No Snakemake imports, so this is directly unit-testable.
"""
import csv
from pathlib import Path

from .phage_inputs import resolve_assemblies, is_compressed


def pgap_dir_to_samplesheet(input_dir, samplesheet_file, taxon='', pattern=None):
    """Scan `input_dir` (recursively) for FASTAs and write name,fasta,taxon."""
    config = {'assembly_dir': str(input_dir), 'recursive': True}
    if pattern:
        config['pattern'] = pattern

    assemblies = resolve_assemblies(config, basedir=Path.cwd())

    compressed = [str(p) for p in assemblies.values() if is_compressed(p)]
    if compressed:
        raise ValueError(
            "PGAP takes FASTA input directly and this workflow has no "
            "decompression step; decompress these first:\n"
            + "\n".join(f"  {p}" for p in sorted(compressed))
        )

    out_dir = Path(samplesheet_file).parent
    if out_dir and not out_dir.exists():
        out_dir.mkdir(parents=True)

    with open(samplesheet_file, 'w', newline='') as fout:
        writer = csv.writer(fout)
        writer.writerow(['name', 'fasta', 'taxon'])
        for name, path in sorted(assemblies.items()):
            writer.writerow([name, str(path), taxon])

    return assemblies
