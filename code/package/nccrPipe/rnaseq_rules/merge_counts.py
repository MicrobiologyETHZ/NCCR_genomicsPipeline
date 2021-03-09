import pandas as pd
import yaml

from pathlib import Path


def parse_config():


    return None

def merge_featureCounts(outdir, samples):
    pd_files = [Path(outdir) / f'{sample}/{sample}.count.txt' for sample in samples]
    pd_list =[]
    for f in pd_files:
        df = pd.read_table(f, skiprows=[0], index_col=0).iloc[:,[4,5]]
        print(df)
        df.columns = ["Length", f.stem.split(".")[0]]
        pd_list.append(df)
    return pd.concat(pd_list, axis=1)


def merge_STARCounts(outdir, samples, strand=0):
    pd_list = [(pd.read_table(Path(outdir)/f'{sample}/{sample}_ReadsPerGene.out.tab', index_col=0, header=None)
                .iloc[:, strand]
                .rename(sample)) for sample in samples]
    return pd.concat(pd_list, axis=1)


def merge_kallistoCounts(outdir, samples):
    pd_list = [(pd.read_table(Path(outdir)/f'{sample}/abundance.tsv',
                             index_col=0)
                .rename({'est_counts': f'{sample}_est_counts',
                         'tpm': f'{sample}_tpm'}, axis=1)) for sample in samples]
    df = pd.concat(pd_list, axis=1)
    return df.loc[:, ~df.columns.duplicated()]

if __name__ == "__main__":
    samples = ['R1', 'R2']
    outdir = 'test'
    print(merge_featureCounts(outdir, samples))