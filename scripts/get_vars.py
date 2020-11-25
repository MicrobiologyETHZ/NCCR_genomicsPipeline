from pathlib import Path
import sys

def get_subsamples(sample_file):
    subsamples = []
    if Path(sample_file).is_file():
        subsamples = set(Path(sample_file).read_text().splitlines())
    if len(subsamples) == 0:
        exit(1)
    return subsamples


def get_samples(dataDir, subsamples, suffix = '.fq.gz'):
    samples = []

    for s in subsamples:
        pattern = f'{s}*{suffix}'
        sample_path = list(Path(dataDir).joinpath(s).rglob(pattern))
        if len(sample_path)<1:
            sys.exit(1)
        else:
            sample = str(sample_path[0].name).split('.')[0]
            samples.append(sample)
    return samples


