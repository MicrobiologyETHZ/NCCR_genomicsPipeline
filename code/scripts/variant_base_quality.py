import sys
from pathlib import Path
import pandas as pd

def get_mpileup_for_SNPs(vcf_file, mpileup_file):
    """
    Take in a vcf file and a mpileup, and subset to only include vcf positions
    :param vcf_file: path to vcf file
    :param mpileup_file: path to mpileup file
    :return: path to the new file
    """
    variants = []
    pileups = []
    with open(vcf_file, 'r') as fh:
        lines = fh.readlines()
        for line in lines:
            if not line.startswith("#"):
                # Get all the genomic locations for the variants
                genome_loc = "\t".join(line.split()[0:2])
                variants.append(genome_loc)

    with open(mpileup_file, 'r') as fh:
        while variants:
            variant = variants.pop(0)
            print(variant)
            line = fh.readline()
            if not line:
                break
            while variant not in line:
                line = fh.readline()
            pileups.append(line)


    with open(Path(mpileup_file).with_suffix(".VCF.mpileup"), 'w') as fo:
        for p in pileups:
            fo.write(p)

    return Path(mpileup_file).with_suffix(".VCF.mpileup")


def read_support(c):
    if c in [',', '.']:
        return 'REF'
    elif c.lower() in ['a', 'c', 'g', 't']:
        return 'ALT'
    else:
        return 'NA'


def parse_pileups(VCF_pileup):
    df_list = []
    with open(VCF_pileup, 'r') as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            data = line.split()
            locID = ["_".join([data[0],data[1]])]*int(data[3])
            print(locID[0])
            snp = list(data[4])
            support = [read_support(c) for c in snp]
            bqO = list(data[5])
            bq = [ord(c)-33 for c in data[5]]
            mq = [ord(c) for c in data[6]]
            bp = [int(i) for i in data[7].split(",")]
            df_list.append(pd.DataFrame([locID, snp, bqO, bq, support, mq, bp]).T)
    pd.concat(df_list).to_csv(Path(VCF_pileup).with_suffix(".csv"))
    return None

if __name__ == "__main__":
    vcf_file = sys.argv[1]
    mpileup_file = sys.argv[2]
    if len(sys.argv) > 3:
        pfile = sys.argv[3]
    else:
        pfile = get_mpileup_for_SNPs(vcf_file, mpileup_file)
    parse_pileups(pfile)