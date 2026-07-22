"""
Merge per-assembly, per-caller phage results into two cross-assembly tables.

Produces:
  phage_predictions.tsv  one row per predicted viral contig, with its caller,
                         length, and CheckV quality/completeness where available
  phage_annotations.tsv  one row per assembly x caller, summarising contig and
                         CDS counts and how many CDSs phynteny rescued from
                         "hypothetical protein"

Every input is treated as optional. Assemblies with no phage produce empty
files upstream (see phage_collect.py), and a prediction-only run has no
annotation outputs at all — both should yield a valid, if sparser, table rather
than an error.

Usage:
    python phage_summary.py --phage-dir OUTDIR/phage --names a,b --callers genomad,cenotetaker \
        --predictions predictions.tsv --annotations annotations.tsv
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

HYPOTHETICAL = 'hypothetical protein'


def read_table(path):
    """Read a TSV, returning an empty frame for missing/empty/headerless files."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep='\t')
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        return pd.DataFrame()


def contig_lengths(fasta):
    """Map contig name -> sequence length, without a Biopython dependency."""
    path = Path(fasta)
    if not path.is_file() or path.stat().st_size == 0:
        return {}

    lengths, contig, length = {}, None, 0
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if contig is not None:
                    lengths[contig] = length
                contig, length = line[1:].strip().split()[0], 0
            else:
                length += len(line.strip())
    if contig is not None:
        lengths[contig] = length
    return lengths


def collect_predictions(phage_dir, names, callers):
    """One row per predicted viral contig, with its position in the genome.

    Joins CheckV quality and the extraction coordinates, so each predicted
    phage can be located on the original assembly rather than only within its
    own extracted sequence.
    """
    rows = []
    for name in names:
        for caller in callers:
            lengths = contig_lengths(
                phage_dir/name/'viral'/f'{name}.{caller}.fna')

            checkv = read_table(
                phage_dir/name/'checkv'/caller/'quality_summary.tsv')
            quality = {}
            if not checkv.empty and 'contig_id' in checkv.columns:
                quality = checkv.set_index('contig_id').to_dict('index')

            coords = read_table(
                phage_dir/name/'viral'/f'{name}.{caller}.coords.tsv')
            located = {}
            if not coords.empty and 'viral_id' in coords.columns:
                located = coords.set_index('viral_id').to_dict('index')

            for contig, length in lengths.items():
                q = quality.get(contig, {})
                loc = located.get(contig, {})
                rows.append({
                    'assembly': name,
                    'caller': caller,
                    'contig': contig,
                    'length': length,
                    'source_contig': loc.get('source_contig'),
                    'source_start': loc.get('source_start'),
                    'source_end': loc.get('source_end'),
                    'is_provirus': loc.get('is_provirus'),
                    'checkv_quality': q.get('checkv_quality'),
                    'completeness': q.get('completeness'),
                    'contamination': q.get('contamination'),
                    'provirus': q.get('provirus'),
                })

    columns = ['assembly', 'caller', 'contig', 'length',
               'source_contig', 'source_start', 'source_end', 'is_provirus',
               'checkv_quality', 'completeness', 'contamination', 'provirus']
    return pd.DataFrame(rows, columns=columns)


def collect_cds(phage_dir, names, callers):
    """Every annotated feature, in original-genome coordinates.

    Concatenates the per-assembly outputs of phage_remap_coords.py.

    Note the `product` column currently carries pharokka's annotation. Enriching
    it with phynteny's refined functions needs the real column names from a
    phynteny run, which have not been observed yet — see the TODO below.
    """
    frames = []
    for name in names:
        for caller in callers:
            table = read_table(
                phage_dir/name/'annotate'/caller
                / f'{name}.{caller}.genome_coords.tsv')
            if not table.empty:
                frames.append(table)

    if not frames:
        return pd.DataFrame(columns=[
            'assembly', 'caller', 'source_contig', 'genome_start', 'genome_end',
            'strand', 'feature_type', 'viral_id', 'extract_start',
            'extract_end', 'feature_id', 'product'])
    return pd.concat(frames, ignore_index=True)


def collect_annotations(phage_dir, names, callers, predictions):
    """One row per assembly x caller, summarising the annotation chain."""
    counts = (predictions.groupby(['assembly', 'caller']).size().to_dict()
              if not predictions.empty else {})

    rows = []
    for name in names:
        for caller in callers:
            cds = read_table(
                phage_dir/name/'annotate'/caller/'phynteny'
                / 'phynteny_per_cds_funcions.tsv')

            n_cds = len(cds)
            n_hypothetical = 0
            if not cds.empty:
                # Column naming varies between phynteny releases; take whichever
                # function-like column is present rather than assuming one.
                func_cols = [c for c in cds.columns
                             if 'function' in c.lower() or 'phrog' in c.lower()]
                if func_cols:
                    n_hypothetical = int(
                        cds[func_cols[0]].astype(str)
                        .str.lower().str.contains(HYPOTHETICAL).sum())

            rows.append({
                'assembly': name,
                'caller': caller,
                'n_viral_contigs': counts.get((name, caller), 0),
                'n_cds': n_cds,
                'n_hypothetical_cds': n_hypothetical,
                'n_assigned_cds': n_cds - n_hypothetical,
            })

    columns = ['assembly', 'caller', 'n_viral_contigs', 'n_cds',
               'n_hypothetical_cds', 'n_assigned_cds']
    return pd.DataFrame(rows, columns=columns)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phage-dir', required=True,
                        help='OUTDIR/phage, containing one dir per assembly.')
    parser.add_argument('--names', required=True,
                        help='Comma-separated assembly names.')
    parser.add_argument('--callers', required=True,
                        help='Comma-separated callers.')
    parser.add_argument('--predictions', required=True,
                        help='Output path for the per-contig table.')
    parser.add_argument('--annotations', required=True,
                        help='Output path for the per-assembly table.')
    parser.add_argument('--cds', required=True,
                        help='Output path for the per-feature genome-coordinate '
                             'table.')
    args = parser.parse_args(argv)

    phage_dir = Path(args.phage_dir)
    names = [n for n in args.names.split(',') if n]
    callers = [c for c in args.callers.split(',') if c]

    predictions = collect_predictions(phage_dir, names, callers)
    annotations = collect_annotations(phage_dir, names, callers, predictions)
    cds = collect_cds(phage_dir, names, callers)

    for path in (args.predictions, args.annotations, args.cds):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(args.predictions, sep='\t', index=False)
    annotations.to_csv(args.annotations, sep='\t', index=False)
    cds.to_csv(args.cds, sep='\t', index=False)

    print(f'Wrote {len(predictions)} predicted contigs -> {args.predictions}')
    print(f'Wrote {len(annotations)} assembly x caller rows -> {args.annotations}')
    print(f'Wrote {len(cds)} features in genome coordinates -> {args.cds}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
