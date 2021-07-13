from pathlib import Path
OUTDIR = Path(config['outDir'])


# todo gzip problems


rule index:
    input: '{file}.{fasta}'
    output: '{file}.{fasta}.bwt',
        marker = touch('{file}.{fasta}.index.done')
    params:
        qerrfile = '{file}.{fasta}.bwa.qerr',
        qoutfile = '{file}.{fasta}.bwa.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/align.yaml'
    shell: "bwa index {input}"



    #todo DO NOT SUBMIT INDEX AND ALIGN in parallel, glitches...




rule align_to_self:
    input:
        fq1 = OUTDIR/'clean_reads/{sample}/{sample}.1.fq.gz',
        fq2 = OUTDIR/'clean_reads/{sample}/{sample}.2.fq.gz',
        fa = OUTDIR/'assembly/{sample}/scaffolds.fasta',
        index = OUTDIR/'assembly/{sample}/scaffolds.fasta.bwt',
        m = OUTDIR/'assembly/{sample}/{sample}.index.done'
    output:
        marker = touch(OUTDIR/'bams/{sample}/{sample}.bwa.done'),
        bam = OUTDIR/'bams/{sample}/{sample}.bam',
        bai = OUTDIR/'bams/{sample}/{sample}.bam.bai'
    params:
        qerrfile = OUTDIR/'bams/{sample}/{sample}.bwa.qerr',
        qoutfile = OUTDIR/'bams/{sample}/{sample}.bwa.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    log:
        log = OUTDIR/'bams/{sample}/{sample}.bwa.log',
        command = OUTDIR/'bams/{sample}/{sample}.bwa.command'
    conda:
        'envs/align.yaml'
    benchmark:
        OUTDIR/'bams/{sample}/{sample}.bwa.benchmark'
    threads:
        8
    shell:
        "bwa mem  -t 4 "
        "-M {input.fa} "
        "{input.fq1} {input.fq2} | samtools sort --reference {input.fa} "
        "-l 9 -@ 4 -O bam -o {output.bam} -T {output.bam}.tmp; samtools index {output.bam}"


if config['reference'].endswith('.fasta') or config['reference'].endswith('.fasta.gz') or config['reference'].endswith('.fna.gz') or config['reference'].endswith('.fna'):
    REFPATH = Path(config['reference']) # todo refactor

elif config['scaffold'] and config['scaffold'].startswith('min'):
    REFPATH = OUTDIR/f'assembly/{config["reference"]}/{config["reference"]}.scaffolds.{config["scaffold"]}.fasta'

    rule gunzip:
        input: str(REFPATH) + ".gz"
        output: REFPATH
        params:
            qerrfile = lambda wildcards: OUTDIR/f'logs/ref_bams/{REFPATH}.gunzip.qerr',
            qoutfile = lambda wildcards: OUTDIR/f'logs/ref_bams/{REFPATH}.gunzip.qout',
            scratch = 6000,
            mem = 7700,
            time = 1400
        log: log= OUTDIR/'logs/ref_bams/gunzip.log'
        threads:
            4
        shell: 'gunzip {input} 2> {log.log}'
else:
    REFPATH = OUTDIR/f'assembly/{config["reference"]}/scaffolds.fasta'

print(REFPATH)

# rule indexFiltered:
#     input: '{sample}.scaffolds.min0.fasta.gz'
#     output: '{sample}.scaffolds.min0.fasta.bwt',
#         marker = touch(OUTDIR/'bams/{sample}/{sample}.indexMin0.done')
#     params:
#         fa = '{sample}.scaffolds.min0.fasta',
#         qerrfile = lambda wildcards: OUTDIR/f'bams/{wildcards.sample}/{wildcards.sample}.indexMin0.qerr',
#         qoutfile = lambda wildcards: OUTDIR/f'bams/{wildcards.sample}/{wildcards.sample}.indexMin0.qout',
#         scratch = 6000,
#         mem = 7700,
#         time = 1400
#     log:
#         log = OUTDIR/'bams/{sample}/{sample}.indexMin0.log',
#         command = OUTDIR/'bams/{sample}/{sample}.indexMin0.command'
#     conda:
#         'envs/align.yaml'
#     benchmark:
#         OUTDIR/'bams/{sample}/{sample}.bwa.benchmark'
#     threads:
#         8
#     shell: 'gunzip {input}; bwa index {params.fa}; gzip {params.fa}'

rule align_to_ref:
    input:
        fq1 = OUTDIR/'clean_reads/{sample}/{sample}.1.fq.gz',
        fq2 = OUTDIR/'clean_reads/{sample}/{sample}.2.fq.gz',
        fa = REFPATH,
        index = str(REFPATH) + '.bwt',
        m = str(REFPATH) + '.index.done'
    output:
        marker = touch(OUTDIR/'ref_bams/{sample}/{sample}.refbwa.done'),
        bam = OUTDIR/'ref_bams/{sample}/{sample}.bam',
        bai = OUTDIR/'ref_bams/{sample}/{sample}.bam.bai'
    params:
        qerrfile = lambda wildcards: OUTDIR/f'ref_bams/{wildcards.sample}/{wildcards.sample}.bwa.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'ref_bams/{wildcards.sample}/{wildcards.sample}.bwa.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    log:
        log = OUTDIR/'ref_bams/{sample}/{sample}.bwa.log',
        command = OUTDIR/'ref_bams/{sample}/{sample}.bwa.command'
    conda:
        'envs/align.yaml'
    benchmark:
        OUTDIR/'ref_bams/{sample}/{sample}.bwa.benchmark'
    threads:
        8
    shell:
        "bwa mem  -t 4 "
        "-M {input.fa} "
        "{input.fq1} {input.fq2} | samtools sort --reference {input.fa} "
        "-l 9 -@ 4 -O bam -o {output.bam} -T {output.bam}.tmp; samtools index {output.bam}"


rule alignment_stats:
    input: OUTDIR/'{bam}/{sample}/{sample}.bam'
    output: OUTDIR/'{bam}/{sample}/{sample}.bam.stats'
    params:
        qerrfile = lambda wildcards: OUTDIR/f'{wildcards.bam}/{wildcards.sample}/{wildcards.sample}.samtools.stats.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'{wildcards.bam}/{wildcards.sample}/{wildcards.sample}.samtools.stats.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    log:
        log = OUTDIR/'{bam}/{sample}/{sample}.samtools.stats.log',
    conda:
        'envs/align.yaml'
    threads:
        8
    shell:
         "samtools stats {input} > {output}"







