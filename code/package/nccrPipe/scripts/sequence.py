import gzip
import Bio.SeqIO.QualityIO
from Bio.Seq import Seq
from Bio import SeqIO
import Bio.SeqIO.FastaIO as FastaIO
from typing import Generator
from pathlib import Path
import sys

#_________________________________
# Parsers
#_________________________________


def snapgeneToFasta(snapfile, fastafile):
    records = SeqIO.parse(snapfile, "snapgene")
    edited_records = []
    for record in records:
        n = 1
        if record.id == '<unknown id>':
            record.id = Path(snapfile).stem + f"_{n}"
        edited_records.append(record)
    count = SeqIO.write(edited_records, fastafile, "fasta")
    print(f'Converted {count} records')
    return fastafile

#___________________________________
class FastA(object):
    '''
    Standard data container for fasta sequences
    '''
    __slots__ = ['header', 'sequence']
    def __init__(self, header: str, sequence: str) -> None:
        self.header = header
        self.sequence = sequence



def revcomp(sequence: str) -> str:
    '''
    Reverse complement a standard nucleotide sequence.
    :param sequence:
    :return:
    '''
    return str(Seq(sequence).reverse_complement())


def stream_fa(sequence_file: str) -> Generator[FastA, None, None]:
    '''
    Read a fastq file either gzipped or not and return it as a stream of tuples
    (Header, Sequence, Quality)
    :param infile:
    :return: Generator[FastA, None, None]
    '''

    if sequence_file.endswith('fq.gz') or sequence_file.endswith('fastq.gz'):
        with gzip.open(sequence_file, 'rt') as handle:
            for header, sequence, qual in Bio.SeqIO.QualityIO.FastqGeneralIterator(handle):
                yield FastA(header, sequence)
    elif sequence_file.endswith('fq') or sequence_file.endswith('fastq'):
        with open(sequence_file) as handle:
            for header, sequence, qual in Bio.SeqIO.QualityIO.FastqGeneralIterator(handle):
                yield FastA(header, sequence)
    elif sequence_file.endswith('fasta.gz') or sequence_file.endswith('fa.gz'):
        with gzip.open(sequence_file, 'rt') as handle:
            for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                yield FastA(header, sequence)
    elif sequence_file.endswith('fasta') or sequence_file.endswith('fa'):
        with open(sequence_file) as handle:
            for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                yield FastA(header, sequence)
    else:
        raise Exception(f'{sequence_file} not a sequence file.')

if __name__ == "__main__":
    files = sys.argv[1:]
    for f in files:
        fastafile = Path(f).parent/(Path(f).stem + '.fasta')
        snapgeneToFasta(f, fastafile)
