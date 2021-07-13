
rule htseq_count:
    input:
        bam = OUTDIR/'ref_bams/{sample}/{sample}.bam',
        gff = config['gffFile'],
    output:
        cnt = OUTDIR/'counts/{sample}/{sample}.tsv',
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.htseqcount.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.htseqcount.qout',
        stranded = 'no',
        feature = 'CDS',
        scratch = 6000,
        mem = 7700,
        time = 1400
    log:
        log = OUTDIR/'logs/{sample}.htseqcount.log',
    conda:
        'envs/count.yaml'
    threads:
        16
    shell:
         'htseq-count -f bam -r pos -m union -s {params.stranded} -t {params.feature} '
         '-i ID {input.bam} {input.gff} --secondary-alignments=ignore '
         '--supplementary-alignments=ignore > {output.cnt}'
