import random

import Bio.SeqIO.QualityIO
from Bio.Seq import Seq
import Bio.SeqIO.FastaIO as FastaIO
from typing import Generator, Tuple, List, Iterator, Dict, DefaultDict, Counter, Set, Union
import sys
#import sequence
import random
import os
from pathlib import Path
import shlex
import subprocess
import sys

def capture(command_str):
    command = shlex.split(command_str)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    return out, err, proc

def to_str(bytes_or_str):
    if isinstance(bytes_or_str, bytes):
        value = bytes_or_str.decode('utf-8')
    else:
        value = bytes_or_str
    return value

def subsample(samples: Tuple, proportions: Tuple, dataDir: str = '/nfs/home/ansintsova/mapStrains/data/raw'):
    assert len(samples) == len(proportions)
    subReadsR1 = '' # todo won't work for large files
    subReadsR2 = ''
    for sample, p in zip(samples, proportions):
        fq1 = Path(list(Path(dataDir).rglob(f'{sample}*.1.fq'))[0])
        fq2 = Path(list(Path(dataDir).rglob(f'{sample}*.2.fq'))[0])
        randomSeed = random.randint(0, 1000)
        command1 = f'seqtk sample -s{randomSeed} {fq1} {p}'
        out1, _, _ = capture(command1)
        command2 = f'seqtk sample -s{randomSeed} {fq2} {p}'
        out2, _, _ = capture(command1)
        subReadsR1 += out1
        subReadsR1 += out2

    return Path('')

subsample(('LL1',), (0.005,))