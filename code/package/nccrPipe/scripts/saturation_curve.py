
import pandas as pd
from pathlib import Path
#import pybedtools
import os
import shlex
import subprocess
import sys

def to_str(bytes_or_str):
    if isinstance(bytes_or_str, bytes):
        value = bytes_or_str.decode("utf-8")
    else:
        value = bytes_or_str
    return value



def to_bytes(bytes_or_str):
    if isinstance(bytes_or_str, str):
        value = bytes_or_str.encode('utf-8')
    else:
        value = bytes_or_str
    return value



# needs bedtools and samtools

class InputError(Exception):
    pass

class Sample:

    def __init__(self, bam_file, gff_file):
        self.bam = Path(bam_file)
        self.gff_file = Path(gff_file)
        self._read_cutoff = 10
        self.subbam = Path('')

    @property
    def cutoff(self):
        return self._read_cutoff

    @cutoff.setter
    def cutoff(self, cutoff):
        self._read_cutoff = cutoff

    def subsample(self, fraction, random_seed=0):

        '''
        :param fraction:
        :param random_seed:
        :return:
        '''
        self.subbam = self.bam.with_suffix(".sub{}.bam".format(str(fraction).split(".")[1]))
        if not 0 < fraction <= 1:
            raise InputError
        cmd = "samtools view -hbs {}.{} {} -t 4 -o {}".format(str(random_seed), str(fraction).split(".")[1],
                                                         self.bam, self.subbam)
        subprocess.run(shlex.split(cmd))
        if self.subbam.exists():
            return self.subbam
        else:
            raise InputError

    def calculate_coverage(self):
        cmd = f"bedtools coverage -a {self.gff_file} -b {self.subbam}"  # todo replace with featureCounts
        output = subprocess.run(shlex.split(cmd), capture_output=True)
        coverage = [s.split("\t") for s in to_str(output.stdout).split("\n") if len(s) > 0]
        df = pd.DataFrame(coverage, columns=["chr", "soft", "feat", "start", "end",
                                            "score", "str", "frame", "attr", "num_reads",
                                            "num_bases_covered", "len_gene", "coverage"])
        return df

    def get_num_expressed_genes(self, reads=True):
        df = self.calculate_coverage()
        if reads:
            df["num_reads"] = df["num_reads"].astype(int)
            return df.loc[df["num_reads"] >= self.cutoff].shape[0]
        else:
            df["coverage"] = df["coverage"].astype(float)
            return df.loc[df["coverage"] >= self.cutoff].shape[0]

    def get_num_expressed_genes_per_fraction(self, fraction, repeat=3, reads=True):
        num_genes = []
        for i in range(repeat):
            print("Subsampling")
            print(i)
            self.subsample(fraction, random_seed=i)
            print('Counting Expressed')
            num_genes.append(self.get_num_expressed_genes(reads))

            os.remove(self.subbam)
        return sum(num_genes)/len(num_genes)

    def get_num_expr_genes_in_fractions(self, fractions, repeat=1, reads=True):
        expr_genes = []
        for fraction in fractions:
            print(fraction)
            fraction_expr = self.get_num_expressed_genes_per_fraction(fraction, repeat, reads)
            expr_genes.append((fraction, fraction_expr))
        return expr_genes


def main():
    bam = sys.argv[1]
    gff = sys.argv[2]
    satFile = Path(bam).with_suffix(".satCurve.csv")
    fractions = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1]
    s = Sample(bam, gff)
    saturation_curve_data = s.get_num_expr_genes_in_fractions(fractions)
    pd.DataFrame(saturation_curve_data, columns=['fraction', 'fraction_expressed']).to_csv(satFile)


if __name__ == "__main__":
    main()

