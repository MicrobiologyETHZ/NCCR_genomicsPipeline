from pathlib import Path

ADAPTERS = Path(config['adapters'])
PHIX = Path(config['phix'])


def getFastq1(wildcards):
    return str(list(Path(DATADIR).joinpath(wildcards.sample).rglob('*.1.fq.gz'))[0])

def getFastq2(wildcards):
    return str(list(Path(DATADIR).joinpath(wildcards.sample).rglob('*.2.fq.gz'))[0])


def return_sample(wildcards):
    return str(Path(OUTDIR/f'{wildcards.subsample}/{wildcards.sample}').parent.stem)



rule qc:
    input:
        fq1 = getFastq1,
        fq2 = getFastq2,
        adapters = ADAPTERS,
        phix = PHIX
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
        trimq=config['trimq'],
        maq=config['mapq'],
        minlen=config['minlen'],
        qoutfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.qc.qout',
        qerrfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.qc.qerr',
        scratch = 500,
        mem = 8000,
        time = 235
    conda:
        "envs/qc.yaml"
    benchmark:
        OUTDIR /'clean_reads/{sample}/{sample}.qc.benchmark'
    log:
        logfile = OUTDIR /'logs/qc/{sample}.qc.log',
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



rule merge_all:
    input:
        fqz1 = OUTDIR /'clean_reads/{sample}/{sample}.1.fq.gz',
        fqz2 = OUTDIR /'clean_reads/{sample}/{sample}.2.fq.gz',
        fqzs = OUTDIR /'clean_reads/{sample}/{sample}.s.fq.gz',
        marker = OUTDIR /'clean_reads/{sample}/{sample}.qc.done'
    output:
        fqz1 = OUTDIR /'merged_reads/{sample}/{sample}.1.fq.gz',
        fqz2 = OUTDIR /'merged_reads/{sample}/{sample}.2.fq.gz',
        fqzs = OUTDIR /'merged_reads/{sample}/{sample}.s.fq.gz',
        fqzm = OUTDIR /'merged_reads/{sample}/{sample}.m.fq.gz',
        hist = OUTDIR /'merged_reads/{sample}/{sample}.merge.hist',
        marker = touch(OUTDIR/'merged_reads/{sample}/{sample}.merge.done')
    params:
        minoverlap = 16,
        scratch = 1000,
        mem = 8000,
        time = 235,
        qerrfile = lambda wildcards: OUTDIR /f'logs/merge/{wildcards.sample}.merge.qerr',
        qoutfile = lambda wildcards: OUTDIR /f'logs/merge/{wildcards.sample}.merge.qout'
    conda:
        "envs/qc.yaml"
    log:
        log = OUTDIR /'logs/merge/{sample}.merge.log',
    benchmark:
        OUTDIR /'merged_reads/{sample}/{sample}.merge.benchmark'
    threads:
        16
    shell:
         "rsync {input.fqzs} {output.fqzs}; "
         "bbmerge.sh -Xmx64G pigz=t bgzip=f threads=4 "
         "overwrite=t in1={input.fqz1} in2={input.fqz2} "
         "out={output.fqzm} outu1={output.fqz1} outu2={output.fqz2} "
         "minoverlap={params.minoverlap} usejni=t ihist={output.hist} &> {log.log}"




## UNDER CONSTRUCTION ##
rule test_fastqc:
    input: '/science/ansintsova/esbl_strains/clean_reads/LL10/LL10_FDSW202435993-1r_HF7HLDSXY_L4.1_fastqc.html'

rule fastqc:
    input: fqz = '{sample}.fq.gz'
    output:
        marker = touch('{sample}.fastqc.done'),
        fqc = '{sample}_fastqc.html'
    threads:
        1
    params:
        outDir = OUTDIR/'QC',
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
        fastqc {input.fqz}  &> {log}
        '''

# todo NOT WORKING right now
# rule multiqc:
#     input:
#         [OUTDIR/f'clean_reads/{sub}/{sam}.1.fastqc.done' for sub, sam in zip(SUBSAMPLES, SAMPLES)][0],
#         [OUTDIR/f'clean_reads/{sub}/{sam}.2.fastqc.done' for sub, sam in zip(SUBSAMPLES, SAMPLES)][0]
#     output:
#         OUTDIR/'QC/multiqc_report.html'
#     params:
#         outdir=OUTDIR/'QC',
#         readDir = OUTDIR/'clean_reads',
#         scratch = 1000,
#         time = 100,
#         mem = 4000,
#         qerrfile = OUTDIR/'QC/QC.multiqc.qerr',
#         qoutfile = OUTDIR/'QC/QC.multiqc.qout'
#     conda:
#         'envs/qc.yaml'
#     shell: "multiqc {params.readDir} -o {params.outdir}"


