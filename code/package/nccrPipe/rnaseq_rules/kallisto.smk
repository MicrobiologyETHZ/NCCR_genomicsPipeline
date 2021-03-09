rule kallisto_index:
    input: tna = config['transcriptome']
    output: marker = touch(OUTDIR/f"{config['ProjectName']}.kallisto.index.done"),
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/kallisto_index.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/kallisto_index.qout',
        kallistoIdx = config["kallistoIdx"],
        scratch = 6000,
        mem = 8000,
        time = 1400
    conda:
        'envs/kallisto.yaml'
    log: OUTDIR/'logs/kallisto_index.log'
    threads:
        32
    shell:
        "kallisto index -i {params.kallistoIdx} {input.tna} "


rule kallisto_quant:
    input: index_done = OUTDIR/f"{config['ProjectName']}.kallisto.index.done",
        fq1 = OUTDIR/'clean_reads/{sample}/{sample}.1.fq.gz',
        fq2= OUTDIR/'clean_reads/{sample}/{sample}.2.fq.gz'
    output: marker = touch(OUTDIR/'kallisto/{sample}.done')
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.kallisto_quant.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.kallisto_quant.qout',
        kallistoIdx = config["kallistoIdx"],
        prefix = lambda wildcards: OUTDIR/f'kallisto/{wildcards.sample}',
        scratch = 6000,
        threads = 8,
        mem = 8000,
        time = 1400
    conda:
        'envs/kallisto.yaml'
    log: OUTDIR/'logs/{sample}.kallisto_quant.log'
    threads:
        32
    shell:
         "kallisto quant -i {params.kallistoIdx} -o {params.prefix} "
         "-t {params.threads} -b 100 <(zcat {input.fq1}) <(zcat {input.fq2}) &> {log}"