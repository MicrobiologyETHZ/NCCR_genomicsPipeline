import Bio.SeqIO.QualityIO
from Bio.Seq import Seq
import Bio.SeqIO.FastaIO as FastaIO
from typing import Generator
import sys
import pandas as pd
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



def GC_skew_per_position(genome: str) -> list:
    skew = [0]
    for base in genome:
        if base == 'G':
            nskew = skew[-1] + 1
        elif base.upper() == 'C':
            nskew = skew[-1] - 1
        else:
            nskew = skew[-1]
        skew.append(nskew)
    return skew

if __name__ == '__main__':
    genome_file = sys.argv[1]
    skew_file = sys.argv[2]
    sequence_generator = stream_fa(genome_file)
    skew = []
    pos = []
    chr = []
    for record in sequence_generator:
        print(record.header)
        print(record.sequence[0:10])
        skew += GC_skew_per_position(record.sequence)
        pos += list(range(len(record.sequence)))
        chr += [record.header.split()[0]]*len(record.sequence)
    fl = [skew[0:10], pos[0:10], chr[0:10]]

    print(fl)
    df = pd.DataFrame(list(zip(skew, pos, chr)), columns = ['skew', 'pos', 'chr'])
    print(df.head())
    df.to_csv(sys.argv[2])





