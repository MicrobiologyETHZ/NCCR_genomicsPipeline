from pathlib import Path
import scripts.get_vars as gv
OUTDIR = Path(config['outDir'])
DATADIR = Path(config['dataDir'])
OUTDIR = Path(config['outDir'])
REFSTRAIN = OUTDIR/'assembly/LL6/LL6.scaffolds.min500.fasta'
SFILE = Path(config['sampleFile'])
SUBSAMPLES = gv.get_subsamples(SFILE)

rule run_quast:
    input: scaf1 = OUTDIR/'{assembly}/{sample}/scaffolds.fasta.gz',
        r1 = OUTDIR/'clean_reads/{sample}/{sample}.1.fq.gz',
        r2 = OUTDIR/'clean_reads/{sample}/{sample}.2.fq.gz',
    output:
        touch(OUTDIR/'{assembly}/{sample}/{sample}.assembly_qc.quast.done'),
        OUTDIR/'{assembly}/{sample}/report.txt'
    params:
        outDir = lambda wildcards: OUTDIR/f'{wildcards.assembly}/{wildcards.sample}',
        qerrfile = OUTDIR/'{assembly}/{sample}/{sample}.quast.qerr',
        qoutfile = OUTDIR/'{assembly}/{sample}/{sample}.quast.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    threads:
        4
    conda:
        'envs/quast.yaml'
    shell:
         'quast.py {input.scaf1} '
         ' -1 {input.r1} -2 {input.r2} '
         '-o {params.outDir} '


rule gzip:
    input: '{sample}.fasta',
    output: '{sample}.fasta.gz',
    params:
        qerrfile = '{sample}.gzip.qerr',
        qoutfile = '{sample}.gzip.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    threads:
        8
    shell:
        'gzip {input}'



# rule run_quast_all:
#     input: [OUTDIR/f'assembly/{sample}/scaffolds.fasta.gz' for sample in SUBSAMPLES]
#     output:
#         touch(OUTDIR/'quast/assembly_qc/assembly_qc.quast.done')
#         #OUTDIR/'quast/report.txt'
#     params:
#         outDir = OUTDIR/'quast/assembly_qc',
#         qerrfile = OUTDIR/'quast/assembly.quast.qerr',
#         qoutfile = OUTDIR/'quast/assemlby.quast.qout',
#         scratch = 6000,
#         mem = 7700,
#         time = 1400
#     threads:
#         4
#     conda:
#         'envs/quast.yaml'
#     shell:
#          'quast.py {input} '
#          '-o {params.outDir} '


