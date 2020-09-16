from pathlib import Path
import scripts.get_vars as gv

DATADIR = Path(config['dataDir'])
OUTDIR = Path(config['outDir'])
REFSTRAIN = OUTDIR/'assembly/LL6/LL6.scaffolds.min500.fasta'
SFILE = Path(config['sampleFile'])
SUBSAMPLES = gv.get_subsamples(SFILE)
SAMPLES = gv.get_samples(DATADIR, SUBSAMPLES)
rule align_genomes:
    input: [OUTDIR/f'mummer/{subsample}/{subsample}.report' for subsample in SUBSAMPLES]

rule dnadiff:
    input:
         ref = REFSTRAIN,
         query = OUTDIR/'assembly/{sample}/{sample}.scaffolds.min500.fasta',
    output:
        OUTDIR/'mummer/{sample}/{sample}.report'
    params:
        sample = '{sample}',
        prefix = OUTDIR/'mummer/',
        qerrfile = OUTDIR/'mummer/{sample}.dnadiff.qerr',
        qoutfile = OUTDIR/'mummer/{sample}.dnadiff.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/compare_genomes.yaml'
    threads:
        8
    shell:
        "dnadiff -p {params.prefix}/{params.sample}/{params.sample} {input.ref} {input.query}"



rule testANI:
    input: OUTDIR/f'ANI/LL6_LL13_fastani.out'

rule calculateANI:
    input: ref = OUTDIR/'assembly/{sample1}/scaffolds.fasta.gz',
         genome = OUTDIR/'assembly/{sample2}/scaffolds.fasta.gz'
    output:
        aniFile = OUTDIR/'ANI/{sample1}_{sample2}_fastani.out'
    params:
        qerrfile = OUTDIR/'mummer/{sample1}_{sample2}.fastani.qerr',
        qoutfile = OUTDIR/'mummer/{sample1}_{sample2}.fastani.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/compare_genomes.yaml'
    threads:
        8
    shell:
        "fastANI -q {input.ref} -r {input.genome} -o {output.aniFile}"





#
# rule gunzip_assembly:
#     input: OUTDIR/'assembly/{sample}.fasta.gz'
#     output: OUTDIR/'assembly/{sample}.fasta'
#     params:
#         qerrfile = OUTDIR/'assembly/{sample}.gzip.qerr',
#         qoutfile = OUTDIR/'assembly/{sample}.gzip.qout',
#         scratch = 6000,
#         mem = 7700,
#         time = 1400
#     conda:
#         'envs/mummer.yaml'
#     threads:
#         8
#     shell:
#         'gunzip {input}'





