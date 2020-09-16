#from pathlib import Path
from scripts.get_vars import *

DATADIR = Path(config['dataDir'])
OUTDIR = Path(config['outDir'])
SAMPLEFILE = Path(config['sampleFile'])
SUBSAMPLES = get_subsamples(SAMPLEFILE)



include: "rules/preprocess.smk"
include: "rules/assemble.smk"
include: "rules/quast.smk"
include: "rules/alignment.smk"
include: "rules/annotate.smk"
include: "rules/compare_genomes.smk"
include: "rules/call_variants.smk"
include: "rules/typing.smk"



rule preprocess:
    input: [OUTDIR/f'clean_reads/{sub}/{sub}.qc.done' for sub in SUBSAMPLES]


rule merge:
    input: [OUTDIR/f'merged_reads/{sam}/{sam}.merge.done' for sam in SUBSAMPLES]


rule assemble:
    input: [OUTDIR/f'assembly/{sub}/{sub}.spades.done' for sub in SUBSAMPLES],
        [OUTDIR/f'plasmid/{sub}/{sub}.spades.done' for sub in SUBSAMPLES]

rule assemble_clean:
    input: [OUTDIR/f'assembly/{sub}/{sub}.assembly_cleanup.done' for sub in SUBSAMPLES],
        [OUTDIR/f'plasmid/{sub}/{sub}.assembly_cleanup.done' for sub in SUBSAMPLES]

rule quast:
    input: [OUTDIR/f'assembly/{sample}/report.txt' for sample in SUBSAMPLES],
        [OUTDIR/f'plasmid/{sample}/report.txt' for sample in SUBSAMPLES]

rule annotate:
    input: [OUTDIR/f'annotation/{sub}/{sub}.prokka.done' for sub in SUBSAMPLES]


rule annotate_plasmid:
    input: [OUTDIR/f'plasmid/annotation/{sub}/{sub}.prokka.done' for sub in SUBSAMPLES]

rule nucmer:
    input: [OUTDIR/f'mummer/{subsample}/{subsample}.report' for subsample in SUBSAMPLES]


rule align:
    input: [OUTDIR/f'bams/{sub}/{sub}.bwa.done' for sub in SUBSAMPLES]

rule align_with_ref:
    input: [OUTDIR/f'ref_bams/{sub}/{sub}.bwa.done' for sub in SUBSAMPLES]


rule call_vars:
    input: [OUTDIR/f'VCF/{sub}/{sub}.vcf.done' for sub in SUBSAMPLES]

rule runANI:
    input: [OUTDIR/f'ANI/{i}_{j}_fastani.out' for i in SUBSAMPLES for j in SUBSAMPLES]

rule markDup:
    input: [OUTDIR/f'bams/{sub}/{sub}.rmdup.bam' for sub in SUBSAMPLES]

# HElPER RULES

rule type:
    input: [OUTDIR/f'typing/{sub}/{sub}.mlst' for sub in SUBSAMPLES]


rule serotype:
    input: [OUTDIR/f'typing/{sub}/output.tsv' for sub in SUBSAMPLES]





# rule gzip:
#     input: '{sample}.fasta',
#     output: '{sample}.fasta.gz',
#     params:
#         qerrfile = '{sample}.gzip.qerr',
#         qoutfile = '{sample}.gzip.qout',
#         scratch = 6000,
#         mem = 7700,
#         time = 1400
#     threads:
#         8
#     shell:
#         'gzip {input}'