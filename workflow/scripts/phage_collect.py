"""
Normalise a caller's viral contigs into a single FASTA.

geNomad and Cenote-Taker 3 write their viral sequences to different places and
with different header conventions. This merges the relevant FASTAs for one
assembly into `<name>.<caller>.fna`, prefixing every contig with the assembly
name so headers stay unique once all assemblies are pooled into the summary
tables.

Missing or empty sources are not an error: an assembly with no phage is a
perfectly normal result, and the downstream rules are written to skip gracefully
when this file is empty.

Usage:
    python phage_collect.py --name SAMPLE --output out.fna src1.fna [src2.fna ...]
"""
import argparse
import sys
from pathlib import Path


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


def collect(name, sources, output):
    """Merge `sources` into `output`, prefixing headers and deduplicating."""
    seen = set()
    written = 0

    with open(output, 'w') as out:
        for source in sources:
            for header, seq in read_fasta(source):
                # geNomad's summary FASTA is documented to already contain
                # proviruses, but the provirus FASTA is passed in as well; drop
                # the repeats rather than depending on which is true.
                contig = header.split()[0] if header else ''
                if contig in seen:
                    continue
                seen.add(contig)
                out.write(f'>{name}|{header}\n')
                out.write('\n'.join(seq) + '\n')
                written += 1

    return written


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--name', required=True,
                        help='Assembly name, used to prefix contig headers.')
    parser.add_argument('--output', required=True, help='Output FASTA path.')
    parser.add_argument('sources', nargs='*',
                        help='Caller FASTAs to merge; missing ones are skipped.')
    args = parser.parse_args(argv)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    written = collect(args.name, args.sources, args.output)

    if written:
        print(f'{args.name}: collected {written} viral contigs -> {args.output}')
    else:
        print(f'{args.name}: no viral contigs found; wrote empty {args.output}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
