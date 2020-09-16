from pathlib import Path
OUTDIR = Path(config['outDir'])


# todo gzip problems


rule index:
    input: '{file}.fasta'
    output: '{file}.fasta.bwt'
    params:
        qerrfile = '{file}.bwa.qerr',
        qoutfile = '{file}.bwa.qout',
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


if config['reference'].endswith('.fasta') or config['reference'].endswith('.fasta.gz'):
    REFPATH = Path(config['reference'])
else:
    REFPATH = OUTDIR/f'assembly/{config["reference"]}/scaffolds.fasta'

rule align_to_ref:
    input:
        fq1 = OUTDIR/'clean_reads/{sample}/{sample}.1.fq.gz',
        fq2 = OUTDIR/'clean_reads/{sample}/{sample}.2.fq.gz',
        fa = REFPATH,
        index = str(REFPATH) + '.bwt'
    output:
        marker = touch(OUTDIR/'ref_bams/{sample}/{sample}.bwa.done'),
        bam = OUTDIR/'ref_bams/{sample}/{sample}.bam',
        bai = OUTDIR/'ref_bams/{sample}/{sample}.bam.bai'
    params:
        qerrfile = OUTDIR/'ref_bams/{sample}/{sample}.bwa.qerr',
        qoutfile = OUTDIR/'ref_bams/{sample}/{sample}.bwa.qout',
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









#rule align_to_ref:
#
# rule gunzip:
#     input: '{sample}.fasta.gz'
#     output: '{sample}.fasta'
#     params:
#         qerrfile = '{sample}.gzip.qerr',
#         qoutfile = '{sample}.gzip.qout',
#         scratch = 6000,
#         mem = 7700,
#         time = 1400
#     threads:
#         8
#     shell:
#         'gunzip {input}'