"""
Normalise a caller's viral contigs into a single FASTA, plus a coordinate map.

Both callers *extract* viral regions from the host contigs, so everything
annotated downstream has coordinates relative to the extract rather than the
genome. This writes the offsets needed to invert that, alongside the sequences.

The two callers use different conventions, verified against a real run:

  geNomad          `_summary/*_virus_summary.tsv`
                   seq_name  = <contig>|provirus_<start>_<end>
                   coordinates = "<start>-<end>", 1-indexed INCLUSIVE
                   (end - start + 1 == length), NA for whole-contig viruses
                   -> offset = start - 1

  Cenote-Taker 3   `*_virus_summary.tsv`  : contig, input_name (original name)
                   `*_prune_summary.tsv`  : contig, chunk_name, chunk_start,
                                            chunk_stop; HALF-OPEN, 0-indexed
                   (chunk_stop - chunk_start == chunk_length)
                   sequences named <contig>@<chunk_name>
                   -> offset = chunk_start

`offset` is normalised to a 0-based genomic start for both, so downstream code
uses one formula: genome_start = offset + extract_start (extract_start 1-based).

Only geNomad's `_summary/*_virus.fna` is used, never `_find_proviruses/`: the
latter holds *candidate* proviral regions, including ones that failed virus
classification. On the reference run it had 9 sequences against the summary's 8.

Missing or empty inputs are not an error — an assembly with no phage is a normal
result, and downstream rules skip gracefully on an empty FASTA.
"""
import argparse
import csv
import sys
from pathlib import Path

COORDS_COLUMNS = [
    'viral_id', 'assembly', 'caller', 'source_contig',
    'offset', 'source_start', 'source_end', 'extracted_length',
    'is_provirus', 'orientation',
]


def read_fasta(path):
    """Yield (header, sequence_lines) pairs. Tolerates a missing file."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return

    header, seq = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    yield header, seq
                header, seq = line[1:], []
            elif header is not None:
                seq.append(line)
    if header is not None:
        yield header, seq


def read_tsv(path):
    """Read a TSV into a list of dicts; empty list if missing or empty."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def genomad_coords(summary_path):
    """Map geNomad seq_name -> extraction coordinates.

    Uses the `coordinates` column, falling back to parsing the name, since both
    encode the same 1-indexed inclusive range.
    """
    coords = {}
    for row in read_tsv(summary_path):
        seq_name = row.get('seq_name', '')
        if not seq_name:
            continue

        raw = (row.get('coordinates') or 'NA').strip()
        source_contig = seq_name.split('|provirus_')[0]

        if raw and raw != 'NA' and '-' in raw:
            start, end = (int(x) for x in raw.split('-', 1))
            is_provirus = True
        else:
            # Whole-contig virus: the extract IS the contig, so no offset.
            start = 1
            end = int(row.get('length') or 0)
            is_provirus = False

        coords[seq_name] = {
            'source_contig': source_contig,
            'offset': start - 1,          # 1-indexed inclusive -> 0-based start
            'source_start': start,
            'source_end': end,
            'is_provirus': is_provirus,
        }
    return coords


def cenotetaker_coords(virus_summary_path, prune_summary_path):
    """Map Cenote-Taker 3 sequence name -> extraction coordinates.

    Needs both files: virus_summary carries the internal-to-original contig
    mapping (`contig` -> `input_name`) but no coordinates; prune_summary carries
    the chunk coordinates but only the internal name.
    """
    # (contig, chunk_name) -> chunk row
    chunks = {}
    for row in read_tsv(prune_summary_path):
        chunks[(row.get('contig', ''), row.get('chunk_name', ''))] = row

    coords = {}
    for row in read_tsv(virus_summary_path):
        seq_name = row.get('contig', '')
        if not seq_name:
            continue

        original = row.get('input_name') or seq_name
        contig, _, chunk_name = seq_name.partition('@')
        length = int(row.get('virus_seq_length') or 0)

        chunk = chunks.get((contig, chunk_name)) if chunk_name else None
        if chunk:
            start = int(chunk['chunk_start'])       # already 0-based
            stop = int(chunk['chunk_stop'])         # half-open
            coords[seq_name] = {
                'source_contig': original,
                'offset': start,
                'source_start': start + 1,          # report 1-based inclusive
                'source_end': stop,
                'is_provirus': True,
            }
        else:
            # Whole contig kept intact — no pruning, so no offset.
            coords[seq_name] = {
                'source_contig': original,
                'offset': 0,
                'source_start': 1,
                'source_end': length,
                'is_provirus': False,
            }
    return coords


def collect(name, caller, fasta_sources, coords, output, coords_out):
    """Write the merged FASTA and its coordinate map.

    Contig headers keep each caller's native ID. They are NOT prefixed with the
    assembly name: geNomad already uses `|` inside provirus names, so prefixing
    produced ambiguous IDs like `Leaf257|contig_1|provirus_1_2`. Uniqueness
    across assemblies comes from the `assembly` column instead.
    """
    seen = set()
    rows = []

    with open(output, 'w') as out:
        for source in fasta_sources:
            for header, seq in read_fasta(source):
                viral_id = header.split()[0] if header else ''
                if not viral_id or viral_id in seen:
                    continue
                seen.add(viral_id)

                out.write(f'>{header}\n')
                out.write('\n'.join(seq) + '\n')

                info = coords.get(viral_id)
                if info is None:
                    # Sequence with no summary row: keep it, but with a zero
                    # offset and a warning, rather than silently dropping it or
                    # inventing coordinates.
                    print(f'WARNING: {viral_id} has no entry in the {caller} '
                          f'summary; genome coordinates will be unavailable',
                          file=sys.stderr)
                    info = {'source_contig': '', 'offset': 0, 'source_start': 0,
                            'source_end': 0, 'is_provirus': False}

                rows.append({
                    'viral_id': viral_id,
                    'assembly': name,
                    'caller': caller,
                    'source_contig': info['source_contig'],
                    'offset': info['offset'],
                    'source_start': info['source_start'],
                    'source_end': info['source_end'],
                    'extracted_length': sum(len(s) for s in seq),
                    'is_provirus': info['is_provirus'],
                    # Both callers were observed to emit forward-oriented
                    # extracts. Carried explicitly so a reverse-complementing
                    # case can be handled without changing the schema.
                    'orientation': '+',
                })

    with open(coords_out, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=COORDS_COLUMNS,
                                delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    return len(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--name', required=True,
                        help='Assembly name.')
    parser.add_argument('--caller', required=True,
                        choices=['genomad', 'cenotetaker'])
    parser.add_argument('--output', required=True, help='Output FASTA path.')
    parser.add_argument('--coords-out', required=True,
                        help='Output coordinate-map TSV path.')
    parser.add_argument('--fasta', required=True,
                        help="The caller's viral sequence FASTA.")
    parser.add_argument('--summary', required=True,
                        help="The caller's virus summary TSV.")
    parser.add_argument('--prune-summary', default=None,
                        help='Cenote-Taker 3 prune summary TSV (required for it).')
    args = parser.parse_args(argv)

    if args.caller == 'genomad':
        coords = genomad_coords(args.summary)
    else:
        if not args.prune_summary:
            parser.error('--prune-summary is required for cenotetaker')
        coords = cenotetaker_coords(args.summary, args.prune_summary)

    for path in (args.output, args.coords_out):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    written = collect(args.name, args.caller, [args.fasta], coords,
                      args.output, args.coords_out)

    if written:
        print(f'{args.name}/{args.caller}: collected {written} viral contigs')
    else:
        print(f'{args.name}/{args.caller}: no viral contigs found; '
              f'wrote empty {args.output}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
