from pathlib import Path
import pandas as pd

DATADIR = Path(config["dataDir"])
OUTDIR = Path(config['outDir'])
import sys

# Assume have a sample sheet:
# Sample Name, Unit, forward reads, reverse reads
sampleInfo = pd.read_csv(config['samples'])

samples_to_merge = (sampleInfo.loc[sampleInfo.groupby('sample')
                    .unit.filter(lambda x: x.nunique() > 1).index]['sample']
                    .unique())



if config['merge_replicates'] and len(samples_to_merge) > 0:

    df1 = sampleInfo[~sampleInfo['sample'].isin(samples_to_merge)][['sample', 'fastq_1', 'fastq_2']]
    df2 = pd.DataFrame([[s, DATADIR/f"concat/{s}/{s}_concat.1.fq.gz", DATADIR/f"concat/{s}/{s}_concat.2.fq.gz"]
                        for s in samples_to_merge], columns=['sample', 'fastq_1', 'fastq_2'])
    sampleInfo_merged = pd.concat([df1, df2])
    sampleInfo_merged['unit'] = ""

else:
    sampleInfo_merged = sampleInfo.copy()

def get_fq1_to_merge(wildcards):
    fq1_to_merge = sampleInfo[sampleInfo['sample'] == wildcards.sample].fastq_1.values
    return fq1_to_merge


def get_fq2_to_merge(wildcards):
    fq2_to_merge = sampleInfo[sampleInfo['sample'] == wildcards.sample].fastq_2.values
    return fq2_to_merge


rule merge_replicates:
    input:
        expand(DATADIR/ "concat/{sample}/{sample}_concat.1.fq.gz", sample=samples_to_merge)


rule concat_fastq:
    input: fq1 =  get_fq1_to_merge,
        fq2 = get_fq2_to_merge
    output: fq1_out = DATADIR/"concat/{sample}/{sample}_concat.1.fq.gz",
        fq2_out = DATADIR/"concat/{sample}/{sample}_concat.2.fq.gz"
    params: fq1 = lambda wildcards: DATADIR/f"concat/{wildcards.sample}/{wildcards.sample}_concat.1.fq",
        fq2 = lambda wildcards: DATADIR/f"concat/{wildcards.sample}/{wildcards.sample}_concat.2.fq",
          # on MacOS use these instead of input
        in1 = lambda wildcards, input: [f"< {f}" for f in input.fq1],
        in2 = lambda wildcards, input: [f"< {f}" for f in input.fq2]
    shell:
        "zcat  {params.in1} > {params.fq1}; gzip {params.fq1}; "
        "zcat  {params.in2} > {params.fq2}; gzip {params.fq2}"


def getFastq1(wildcards):
    return sampleInfo_merged[sampleInfo_merged['sample'] == wildcards.sample].fastq_1

def getFastq2(wildcards):
    return sampleInfo_merged[sampleInfo_merged['sample'] == wildcards.sample].fastq_2


if config['qc']:
    if not config['se']:
        rule qc:
            input:
                fq1 = getFastq1,
                fq2 = getFastq2,
                adapters = Path(config['adapters']),
                phix = Path(config['phix'])
            output:
                fq1_clean = OUTDIR /'clean_reads/{sample}/{sample}.1.fq.gz',
                fq2_clean = OUTDIR /'clean_reads/{sample}/{sample}.2.fq.gz',
                adapter_matched = OUTDIR /'clean_reads/{sample}/removedreads/{sample}.adapter.matched.fq.gz',
                adapter_singletons = OUTDIR /'clean_reads/{sample}/removedreads/{sample}.adapter.singletons.fq.gz',
                adapter_stats = OUTDIR /'clean_reads/{sample}/{sample}.adapter.stats',
                phix_matched = OUTDIR /'clean_reads/{sample}/removedreads/{sample}.phix.matched.fq.gz',
                phix_singletons = OUTDIR /'clean_reads/{sample}/removedreads/{sample}.phix.singletons.fq.gz',
                phix_stats = OUTDIR /'clean_reads/{sample}/{sample}.phix.stats',
                qc_failed = OUTDIR /'clean_reads/{sample}/removedreads/{sample}.qc.failed.fq.gz',
                qc_singletons = OUTDIR /'clean_reads/{sample}/{sample}.s.fq.gz',
                qc_stats = OUTDIR /'clean_reads/{sample}/{sample}.qc.stats',
                marker = touch(OUTDIR /'clean_reads/{sample}/{sample}.qc.done')
            params:
                trimq = config['trimq'],
                maq = config['mapq'],
                minlen = config['minlen'],
                qoutfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.qc.qout',
                qerrfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.qc.qerr',
                scratch = 500,
                mem = 8000,
                time = 235
            conda:
                "preprocessing"
            benchmark:
                OUTDIR /'clean_reads/{sample}/{sample}.qc.benchmark'
            log:
                log = OUTDIR /'logs/qc/{sample}.qc.log'
            threads:
                8
            shell:
                "bbduk.sh -Xmx1G pigz=t bgzip=f usejni=t "
                "in={input.fq1} in2={input.fq2} "
                "out=stdout.fq outm={output.adapter_matched} "
                "outs={output.adapter_singletons} "
                "refstats={output.adapter_stats} statscolumns=5 "
                "overwrite=t ref={input.adapters} "
                "ktrim=r k=23 mink=11 hdist=1  2>> {log.log} | "
                "bbduk.sh -Xmx1G usejni=t pigz=t bgzip=f "
                "interleaved=true overwrite=t "
                "in=stdin.fq out=stdout.fq "
                "outm={output.phix_matched} outs={output.phix_singletons} "
                "ref={input.phix} k=31 hdist=1 "
                "refstats={output.phix_stats} statscolumns=5 2>> {log.log}| "
                "bbduk.sh -Xmx1G pigz=t bgzip=f usejni=t  "
                "overwrite=t interleaved=true "
                "in=stdin.fq fastawrap=10000 "
                "out1={output.fq1_clean} out2={output.fq2_clean} "
                "outm={output.qc_failed} outs={output.qc_singletons} "
                "minlength={params.minlen} qtrim=rl maq={params.maq} maxns=1  "
                "stats={output.qc_stats} statscolumns=5 "
                "trimq={params.trimq}  2>> {log.log};"
    else:
        rule qc:
            input:
                fq1=getFastq1,
                adapters=Path(config['adapters']),
                phix=Path(config['phix'])
            output:
                fq1_clean=OUTDIR / 'clean_reads/{sample}/{sample}.1.fq.gz',
                adapter_matched=OUTDIR / 'clean_reads/{sample}/removedreads/{sample}.adapter.matched.fq.gz',
                adapter_stats=OUTDIR / 'clean_reads/{sample}/{sample}.adapter.stats',
                phix_matched=OUTDIR / 'clean_reads/{sample}/removedreads/{sample}.phix.matched.fq.gz',
                phix_stats=OUTDIR / 'clean_reads/{sample}/{sample}.phix.stats',
                qc_failed=OUTDIR / 'clean_reads/{sample}/removedreads/{sample}.qc.failed.fq.gz',
                qc_stats=OUTDIR / 'clean_reads/{sample}/{sample}.qc.stats',
                marker=touch(OUTDIR / 'clean_reads/{sample}/{sample}.qc.done')
            params:
                trimq=config['trimq'],
                maq=config['mapq'],
                minlen=config['minlen'],
                qoutfile=lambda wildcards: OUTDIR / f'logs/qc/{wildcards.sample}.qc.qout',
                qerrfile=lambda wildcards: OUTDIR / f'logs/qc/{wildcards.sample}.qc.qerr',
                scratch=500,
                mem=8000,
                time=235
            conda:
                "preprocessing"
            log:
                log=OUTDIR / 'logs/{sample}.qc.log'
            threads:
                8
            shell:
                "bbduk.sh -Xmx1G pigz=t bgzip=f usejni=t "
                "in={input.fq1} out=stdout.fq outm={output.adapter_matched} "
                "refstats={output.adapter_stats} statscolumns=5 "
                "overwrite=t ref={input.adapters} "
                "ktrim=r k=23 mink=11 hdist=1  2>> {log.log} | "
                "bbduk.sh -Xmx1G usejni=t pigz=t bgzip=f "
                "overwrite=t interleaved=f in=stdin.fq out=stdout.fq "
                "outm={output.phix_matched} ref={input.phix} k=31 hdist=1 "
                "refstats={output.phix_stats} statscolumns=5 2>> {log.log}| "
                "bbduk.sh -Xmx1G pigz=t bgzip=f usejni=t  "
                "overwrite=t interleaved=f in=stdin.fq fastawrap=10000 "
                "out={output.fq1_clean} outm={output.qc_failed} "
                "minlength={params.minlen} qtrim=rl maq={params.maq} maxns=1  "
                "stats={output.qc_stats} statscolumns=5 "
                "trimq={params.trimq}  2>> {log.log};"

else:
    rule qc:
        input:
            fq1 = getFastq1,
            fq2 = getFastq2,
        output:
            fq1_clean = OUTDIR /'clean_reads/{sample}/{sample}.1.fq.gz',
            fq2_clean = OUTDIR /'clean_reads/{sample}/{sample}.2.fq.gz',
            marker = touch(OUTDIR /'clean_reads/{sample}/{sample}.qc.done')
        params:
            qoutfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.qc.qout',
            qerrfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.qc.qerr',
            scratch = 500,
            mem = 8000,
            time = 235
        conda:
            "preprocessing"
        benchmark:
            OUTDIR /'clean_reads/{sample}/{sample}.qc.benchmark'
        log:
            log = OUTDIR /'logs/qc/{sample}.qc.log'
        threads:
            8
        shell:
            "cp {input.fq1} {output.fq1_clean}; cp {input.fq2} {output.fq2_clean} "



rule fastqc_clean:
        input:
            fq1 = OUTDIR / 'clean_reads/{sample}/{sample}.1.fq.gz',
            fq2 = OUTDIR / 'clean_reads/{sample}/{sample}.2.fq.gz'
        output:
            fasqc_marker = touch(OUTDIR /'fastqc_clean/{sample}.fastqc_clean.done'),
        params:
            outdir = OUTDIR/'fastqc_clean',
            qoutfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.fastqc_clean.qout',
            qerrfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.fastqc_clean.qerr',
            scratch = 500,
            mem = 8000,
            time = 235
        conda:
            "preprocessing"
        log:
            log = OUTDIR /'logs/qc/{sample}.fastqc_clean.log'
        threads:
            8
        shell:
            "fastqc -o {params.outdir} --noextract {input.fq1} {input.fq2} "


rule fastqc:
        input:
            fq1 = getFastq1,
            fq2 = getFastq2,
        output:
            fasqc_marker = touch(OUTDIR /'fastqc/{sample}.fastqc.done'),
        params:
            outdir = OUTDIR/'fastqc',
            qoutfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.fastqc.qout',
            qerrfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.fastqc.qerr',
            scratch = 500,
            mem = 8000,
            time = 235
        conda:
            "preprocessing"
        log:
            log = OUTDIR /'logs/qc/{sample}.fastqc.log'
        threads:
            8
        shell:
            "fastqc -o {params.outdir} --noextract {input.fq1} {input.fq2} "