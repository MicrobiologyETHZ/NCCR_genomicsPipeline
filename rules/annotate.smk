from pathlib import Path
import scripts.get_vars as gv

DATADIR = Path(config['dataDir'])
OUTDIR = Path(config['outDir'])
ADAPTERS = Path(config['adapters'])
PHIX = Path(config['phix'])
SFILE = Path(config['sampleFile'])

SUBSAMPLES = gv.get_subsamples(SFILE)
SAMPLES = gv.get_samples(DATADIR, SUBSAMPLES)


rule annotate:
    input: [OUTDIR/f'annotation/{sub}/{sam}.prokka.done' for sub, sam in zip(SUBSAMPLES, SAMPLES)][0]


rule prokka:
    input:
        scaffolds = OUTDIR/'assembly/{subsample}/scaffolds.fasta',
        marker = OUTDIR/'assembly/{subsample}/{sample}.spades.done'
    output:
        gff = OUTDIR/'annotation/{subsample}/{sample}.gff',
        gbk = OUTDIR/'annotation/{subsample}/{sample}.gbk',
        fna = OUTDIR/'annotation/{subsample}/{sample}.fna',
        faa = OUTDIR/'annotation/{subsample}/{sample}.faa',
        marker = touch(OUTDIR/'annotation/{subsample}/{sample}.prokka.done')
    params:
        locustag = '{subsample}',
        outdir = OUTDIR/'annotation',
        scratch = 1000,
        mem = 4000,
        time = 235,
        qerrfile = OUTDIR/'annotation/{subsample}/{sample}.prokka.qerr',
        qoutfile = OUTDIR/'annotation/{subsample}/{sample}.prokka.qout'
    conda:
        'envs/annotation.yaml'
    threads:
        2
    shell:
        'prokka --outdir {params.outdir}/{params.locustag} '
        '--locustag {params.locustag} '
        '--prefix {params.locustag} {input.scaffolds} '
#'--proteins NewToxins.faa   '


# rule prodigal:
#     input:
#         scaffolds = '{path}/'+SPADES_FOLDER_NAME+'/{sample}.scaffolds.min500.fasta.gz',
#         marker = '{path}/'+SPADES_FOLDER_NAME+'/{sample}.assembly_cleanup.done'
#     output:
#         faagz = '{path}/'+PRODIGAL_FOLDER_NAME+'/{sample}.scaffolds.min500.genes.faa.gz',
#         fnagz = '{path}/'+PRODIGAL_FOLDER_NAME+'/{sample}.scaffolds.min500.genes.fna.gz',
#         gffgz = '{path}/'+PRODIGAL_FOLDER_NAME+'/{sample}.scaffolds.min500.genes.gff.gz',
#         marker = touch('{path}/'+PRODIGAL_FOLDER_NAME+'/{sample}.scaffolds.min500.prodigal.done')
#     params:
#         faa = '{path}/'+PRODIGAL_FOLDER_NAME+'/{sample}.scaffolds.min500.genes.faa',
#         fna = '{path}/'+PRODIGAL_FOLDER_NAME+'/{sample}.scaffolds.min500.genes.fna',
#         gff = '{path}/'+PRODIGAL_FOLDER_NAME+'/{sample}.scaffolds.min500.genes.gff',
#         scratch = 1000,
#         mem = 4000,
#         time = 235,
#         qerrfile = '{path}/'+PRODIGAL_FOLDER_NAME+'/{sample}.scaffolds.min500.prodigal.qerr',
#         qoutfile = '{path}/'+PRODIGAL_FOLDER_NAME+'/{sample}.scaffolds.min500.prodigal.qout'
#     threads:
#         2
#     log:
#         command = '{path}/'+PRODIGAL_FOLDER_NAME+'/{sample}.scaffolds.min500.prodigal.command',
#     shell:
#         '''
#         #!/bin/bash
#         command="
#         zcat {input.scaffolds} | prodigal -a {params.faa} -d {params.fna} -f gff -o {params.gff} -c -q -m -p meta;
#         pigz -p {threads} {params.faa};
#         pigz -p {threads} {params.fna};
#         pigz -p {threads} {params.gff};
#         ";
#         echo "$command" > {log.command};
#         eval "$command"
#         '''

