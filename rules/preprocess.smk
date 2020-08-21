from pathlib import Path
import sys
DATADIR = Path(config['dataDir'])
OUTDIR = Path(config['outDir'])
ADAPTERS = Path(config['adapters'])
PHIX = Path(config['phix'])
SFILE = Path(config['sampleFile'])

def get_subsamples(sample_file=SFILE):
    subsamples = []
    if Path(sample_file).is_file():
        subsamples = set(Path(sample_file).read_text().splitlines())
    if len(subsamples) == 0:
        exit(1)
    return subsamples

def get_samples(subsamples, suffix = '.fq.gz'):
    samples = []

    for s in subsamples:
        pattern = f'{s}*{suffix}'
        sample_path = list(DATADIR.joinpath(s).rglob(pattern))
        if len(sample_path)<1:
            sys.exit(1)
        else:
            sample = str(sample_path[0].name).split('.')[0]
            samples.append(sample)
    return samples


# def getFastqPair(wildcards):
#     fq1 = Path(DATADIR).joinpath(wildcards.sample).rglob('*.1.fq.gz')
#     fq2 = Path(DATADIR).joinpath(wildcards.sample).rglob('*.2.fq.gz')
#     return fq1, fq2


SUBSAMPLES = get_subsamples()
SAMPLES = get_samples(SUBSAMPLES)

# todo change PE rules
rule qc_samples:
    input: [OUTDIR/f'merged_reads/{sub}/{sam}.merge.done' for sub, sam in zip(SUBSAMPLES, SAMPLES)]


rule qc:
    input:
        fq1 = DATADIR/'{subsample}/{sample}.1.fq.gz',
        fq2 = DATADIR/'{subsample}/{sample}.2.fq.gz',
        adapters = ADAPTERS,
        phix = PHIX
    output:
        fq1_clean = OUTDIR /'clean_reads/{subsample}/{sample}.1.fq.gz',
        fq2_clean = OUTDIR /'clean_reads/{subsample}/{sample}.2.fq.gz',
        adapter_matched = OUTDIR /'clean_reads/{subsample}/removedreads/{sample}.adapter.matched.fq.gz',
        adapter_singletons = OUTDIR /'clean_reads/{subsample}/removedreads/{sample}.adapter.singletons.fq.gz',
        adapter_stats = OUTDIR /'clean_reads/{subsample}/{sample}.adapter.stats',
        phix_matched = OUTDIR /'clean_reads/{subsample}/removedreads/{sample}.phix.matched.fq.gz',
        phix_singletons = OUTDIR /'clean_reads/{subsample}/removedreads/{sample}.phix.singletons.fq.gz',
        phix_stats = OUTDIR /'clean_reads/{subsample}/{sample}.phix.stats',
        qc_failed = OUTDIR /'clean_reads/{subsample}/removedreads/{sample}.qc.failed.fq.gz',
        qc_singletons = OUTDIR /'clean_reads/{subsample}/{sample}.s.fq.gz',
        qc_stats = OUTDIR /'clean_reads/{subsample}/{sample}.qc.stats',
        marker = touch(OUTDIR /'clean_reads/{subsample}/{sample}.qc.done')
    params:
        trimq=config['trimq'],
        maq=config['mapq'],
        minlen=config['minlen'],
        qoutfile = OUTDIR /'clean_reads/{subsample}/{sample}.qc.qout',
        qerrfile = OUTDIR /'clean_reads/{subsample}/{sample}.qc.qerr',
        scratch = 500,
        mem = 8000,
        time = 235

    conda:
        "envs/qc.yaml"
    benchmark:
        OUTDIR /'clean_reads/{subsample}/{sample}.qc.benchmark'
    log:
        logfile = OUTDIR /'clean_reads/{subsample}/{sample}.qc.log',
        command = OUTDIR /'clean_reads/{subsample}/{sample}.qc.command'
    threads:
        8
    shell:
        "bbduk.sh -Xmx1G pigz=t bgzip=f usejni=t " # not sure what this does
        "in={input.fq1} in2={input.fq2} "
        "out=stdout.fq outm={output.adapter_matched} "
        "outs={output.adapter_singletons} "
        "refstats={output.adapter_stats} statscolumns=5 "
        "overwrite=t ref={input.adapters} "
        "ktrim=r k=23 mink=11 hdist=1  2>> {log.logfile} | "
        "bbduk.sh -Xmx1G usejni=t pigz=t bgzip=f "
        "interleaved=true overwrite=t "
        "in=stdin.fq out=stdout.fq "
        "outm={output.phix_matched} outs={output.phix_singletons} "
        "ref={input.phix} k=31 hdist=1 "
        "refstats={output.phix_stats} statscolumns=5 2>> {log.logfile}| "
        "bbduk.sh -Xmx1G pigz=t bgzip=f usejni=t  "
        "overwrite=t interleaved=true "
        "in=stdin.fq fastawrap=10000 "
        "out1={output.fq1_clean} out2={output.fq2_clean} "
        "outm={output.qc_failed} outs={output.qc_singletons} "
        "minlength={params.minlen} qtrim=rl maq={params.maq} maxns=1  "
        "stats={output.qc_stats} statscolumns=5 "
        "trimq={params.trimq}  2>> {log.logfile};"



rule merge:
    input:
        fqz1 = OUTDIR /'clean_reads/{subsample}/{sample}.1.fq.gz',
        fqz2 = OUTDIR /'clean_reads/{subsample}/{sample}.2.fq.gz',
        fqzs = OUTDIR /'clean_reads/{subsample}/{sample}.s.fq.gz',
        marker = OUTDIR /'clean_reads/{subsample}/{sample}.qc.done'
    output:
        fqz1 = OUTDIR /'merged_reads/{subsample}/{sample}.1.fq.gz',
        fqz2 = OUTDIR /'merged_reads/{subsample}/{sample}.2.fq.gz',
        fqzs = OUTDIR /'merged_reads/{subsample}/{sample}.s.fq.gz',
        fqzm = OUTDIR /'merged_reads/{subsample}/{sample}.m.fq.gz',
        hist = OUTDIR /'merged_reads/{subsample}/{sample}.merge.hist',
        marker = touch(OUTDIR/'merged_reads/{subsample}/{sample}.merge.done')
    params:
        minoverlap = 16,
        scratch = 1000,
        mem = 8000,
        time = 235,
        qerrfile = OUTDIR /'merged_reads/{subsample}/{sample}.merge.qerr',
        qoutfile = OUTDIR /'merged_reads/{subsample}/{sample}.merge.qout'
    conda:
        "envs/qc.yaml"
    log:
        log = OUTDIR /'merged_reads/{subsample}/{sample}.merge.log',
        command = OUTDIR /'merged_reads/{subsample}/{sample}.merge.command'
    benchmark:
        OUTDIR /'merged_reads/{subsample}/{sample}.merge.benchmark'
    threads:
        16
    shell:
         "rsync {input.fqzs} {output.fqzs}; "
         "bbmerge.sh -Xmx64G pigz=t bgzip=f threads=4 "
         "overwrite=t in1={input.fqz1} in2={input.fqz2} "
         "out={output.fqzm} outu1={output.fqz1} outu2={output.fqz2} "
         "minoverlap={params.minoverlap} usejni=t ihist={output.hist} &> {log.log}"



rule fastqc:
    input: fqz = '{sample}.fq.gz'
    output:
        marker = touch('{sample}.fastqc.done'),
        fqc = '{sample}_fastqc.html'
    threads:
        1
    params:
        scratch = 1000,
        time = 100,
        mem = 4000,
        qerrfile = '{sample}.fastqc.qerr',
        qoutfile = '{sample}.fastqc.qout'
    conda:
        'envs/qc.yaml'
    log:
        '{sample}.fastqc.log'
    shell:
        '''
        fastqc {input.fqz} &> {log}
        '''


rule multiqc:
    input:
        [OUTDIR/f'clean_reads/{sub}/{sam}.1.fastqc.done' for sub, sam in zip(SUBSAMPLES, SAMPLES)],
        [OUTDIR/f'clean_reads/{sub}/{sam}.2.fastqc.done' for sub, sam in zip(SUBSAMPLES, SAMPLES)]
    output:
        OUTDIR/'QC/multiqc_report.html'
    params:
        outdir=OUTDIR/'QC',
        scratch = 1000,
        time = 100,
        mem = 4000,
    conda:
        'envs/qc.yaml'
    shell: "multiqc {params.outdir} -o {params.outdir}"



        # '''
        #  #!/bin/bash
        #  command="
        #  rsync {input.fqzs} {output.fqzs}
        #  bbmerge.sh -Xmx64G pigz=t bgzip=f threads={threads} overwrite=t in1={input.fqz1} in2={input.fqz2} out={output.fqzm} outu1={output.fqz1} outu2={output.fqz2} minoverlap={params.minoverlap} usejni=t ihist={output.hist} &> {log.log}
        #  ";
        #  echo "$command" > {log.command};
        #  eval "$command"
        # '''