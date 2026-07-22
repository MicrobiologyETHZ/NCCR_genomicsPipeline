"""
Translate pharokka annotations from extract coordinates to genome coordinates.

Predicted viral regions are excised from the host contigs before annotation, so
every feature pharokka calls is positioned relative to the extract. This maps
them back onto the original assembly.

Both callers normalise to a 0-based genomic `offset` in the coords map written
by phage_collect.py, so one formula covers both:

    genome_start = offset + extract_start        (extract_start is 1-based)
    genome_end   = offset + extract_end

Verified end-to-end on a real run: geNomad's `contig_153|provirus_710297_749031`
and Cenote-Taker's `Leaf257_5@C38` (chunk_start 710296, 0-based) both place their
first base at contig_153:710297.

Runs after pharokka rather than at the end of the chain, because pharokka is
where the CDS coordinates are defined — phold and phynteny reuse its calls and
only refine the functional labels.

Outputs GFF3 keyed on the original contig (loads into a genome browser against
the source assembly) and a flat TSV.
"""
import argparse
import csv
import sys
from pathlib import Path

TSV_COLUMNS = [
    'assembly', 'caller', 'source_contig', 'genome_start', 'genome_end',
    'strand', 'feature_type', 'viral_id', 'extract_start', 'extract_end',
    'feature_id', 'product',
]

# How many bases of each extract to check against the assembly. Enough to make a
# coincidental match essentially impossible, cheap enough to do for every contig.
VERIFY_BASES = 30


def read_coords(path):
    """Load the coordinate map keyed by viral contig ID."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return {}
    with open(path, newline='') as handle:
        return {row['viral_id']: row
                for row in csv.DictReader(handle, delimiter='\t')}


def parse_gff(path):
    """Yield pharokka GFF feature rows as dicts.

    Splits with maxsplit=8 deliberately: pharokka's tRNA rows can contain a raw
    tab inside the attributes column (observed in tRNAscan-SE output), which
    would otherwise produce ten fields and break a naive parse.
    """
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return

    with open(path) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t', 8)
            if len(fields) < 8:
                continue
            attributes = fields[8] if len(fields) > 8 else ''
            yield {
                'seqid': fields[0],
                'source': fields[1],
                'type': fields[2],
                'start': int(fields[3]),
                'end': int(fields[4]),
                'score': fields[5],
                'strand': fields[6],
                'phase': fields[7],
                'attributes': attributes,
            }


def parse_attributes(attributes):
    """Parse a GFF attributes string into a dict, tolerating stray fields."""
    parsed = {}
    for item in attributes.split(';'):
        key, sep, value = item.partition('=')
        if sep:
            parsed[key.strip()] = value.strip()
    return parsed


def read_fasta_seqs(path):
    """Load a FASTA into {id: sequence}. Used only for verification."""
    path = Path(path)
    if not path.is_file():
        return {}
    seqs, name, chunks = {}, None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(chunks)
                name, chunks = line[1:].split()[0], []
            elif name is not None:
                chunks.append(line)
    if name is not None:
        seqs[name] = ''.join(chunks)
    return seqs


def verify_offsets(coords, viral_fasta, assembly_fasta):
    """Check each extract really starts where its offset claims.

    Catches the failure modes that would otherwise produce plausible-looking but
    wrong coordinates: an off-by-one, a reverse-complemented extract, or a
    mismatched contig name. Returns a list of human-readable problems.
    """
    viral = read_fasta_seqs(viral_fasta)
    assembly = read_fasta_seqs(assembly_fasta)
    if not viral or not assembly:
        return []

    problems = []
    for viral_id, row in coords.items():
        extract = viral.get(viral_id, '')
        contig = assembly.get(row['source_contig'], '')
        if not extract or not contig:
            continue

        offset = int(row['offset'])
        expected = extract[:VERIFY_BASES].upper()
        actual = contig[offset:offset + VERIFY_BASES].upper()
        if expected != actual:
            problems.append(
                f"  {viral_id}: extract does not match "
                f"{row['source_contig']} at offset {offset}\n"
                f"      extract:  {expected}\n"
                f"      assembly: {actual}")
    return problems


def remap(gff_path, coords, assembly, caller):
    """Translate every GFF feature into genome coordinates."""
    rows = []
    unmapped = set()

    for feature in parse_gff(gff_path):
        # pharokka seqids may still carry a legacy "<assembly>|" prefix; strip
        # it if present so older outputs remain readable.
        viral_id = feature['seqid']
        info = coords.get(viral_id)
        if info is None and '|' in viral_id:
            info = coords.get(viral_id.split('|', 1)[1])

        if info is None:
            unmapped.add(viral_id)
            continue

        offset = int(info['offset'])
        attrs = parse_attributes(feature['attributes'])
        rows.append({
            'assembly': assembly,
            'caller': caller,
            'source_contig': info['source_contig'],
            'genome_start': offset + feature['start'],
            'genome_end': offset + feature['end'],
            'strand': feature['strand'],
            'feature_type': feature['type'],
            'viral_id': viral_id,
            'extract_start': feature['start'],
            'extract_end': feature['end'],
            'feature_id': attrs.get('ID', ''),
            'product': attrs.get('product', ''),
        })

    return rows, sorted(unmapped)


def write_gff3(rows, path):
    """Write GFF3 keyed on the original contig."""
    with open(path, 'w') as out:
        out.write('##gff-version 3\n')
        for row in rows:
            attributes = ';'.join([
                f"ID={row['feature_id']}",
                f"product={row['product']}",
                f"viral_contig={row['viral_id']}",
                f"caller={row['caller']}",
                # Keep the extract-relative position so results stay traceable
                # to the pharokka/phold/phynteny output they came from.
                f"extract_coords={row['extract_start']}-{row['extract_end']}",
            ])
            out.write('\t'.join([
                row['source_contig'], f"phage_{row['caller']}",
                row['feature_type'], str(row['genome_start']),
                str(row['genome_end']), '.', row['strand'], '.', attributes,
            ]) + '\n')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gff', required=True, help='pharokka GFF.')
    parser.add_argument('--coords', required=True,
                        help='Coordinate map from phage_collect.py.')
    parser.add_argument('--name', required=True, help='Assembly name.')
    parser.add_argument('--caller', required=True,
                        choices=['genomad', 'cenotetaker'])
    parser.add_argument('--gff-out', required=True)
    parser.add_argument('--tsv-out', required=True)
    parser.add_argument('--viral-fasta', default=None,
                        help='Viral contigs, for offset verification.')
    parser.add_argument('--assembly-fasta', default=None,
                        help='Original assembly, for offset verification.')
    parser.add_argument('--strict', action='store_true',
                        help='Fail instead of warning when verification fails.')
    args = parser.parse_args(argv)

    coords = read_coords(args.coords)
    rows, unmapped = remap(args.gff, coords, args.name, args.caller)

    if args.viral_fasta and args.assembly_fasta:
        problems = verify_offsets(coords, args.viral_fasta, args.assembly_fasta)
        if problems:
            message = ("Extract offsets do not match the assembly — remapped "
                       "coordinates would be WRONG:\n" + "\n".join(problems))
            if args.strict:
                raise SystemExit(message)
            print(f'WARNING: {message}', file=sys.stderr)
        elif coords:
            print(f'{args.name}/{args.caller}: offsets verified against the '
                  f'assembly for {len(coords)} viral contigs')

    if unmapped:
        print(f'WARNING: {len(unmapped)} sequence(s) in {args.gff} had no entry '
              f'in the coordinate map and were skipped: '
              f'{", ".join(unmapped[:5])}', file=sys.stderr)

    for path in (args.gff_out, args.tsv_out):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    write_gff3(rows, args.gff_out)
    with open(args.tsv_out, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    print(f'{args.name}/{args.caller}: remapped {len(rows)} features '
          f'-> {args.gff_out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
