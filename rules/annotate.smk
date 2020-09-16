from pathlib import Path
import scripts.get_vars as gv

DATADIR = Path(config['dataDir'])
OUTDIR = Path(config['outDir'])
ADAPTERS = Path(config['adapters'])
PHIX = Path(config['phix'])
SFILE = Path(config['sampleFile'])

SUBSAMPLES = gv.get_subsamples(SFILE)
SAMPLES = gv.get_samples(DATADIR, SUBSAMPLES)


rule annotate_all:
    input: [OUTDIR/f'annotation/{sub}/{sub}.prokka.done' for sub in SUBSAMPLES][0],

rule zip:
    input: [OUTDIR/f'assembly/{sample}/scaffolds.fasta.gz' for sample in SUBSAMPLES]


rule prokka:
    input:
        scaffolds = OUTDIR/'assembly/{sample}/scaffolds.fasta',
        marker = OUTDIR/'assembly/{sample}/{sample}.spades.done'
    output:
        gff = OUTDIR/'annotation/{sample}/{sample}.gff',
        gbk = OUTDIR/'annotation/{sample}/{sample}.gbk',
        fna = OUTDIR/'annotation/{sample}/{sample}.fna',
        faa = OUTDIR/'annotation/{sample}/{sample}.faa',
        marker = touch(OUTDIR/'annotation/{sample}/{sample}.prokka.done')
    params:
        locustag = '{sample}',
        outdir = OUTDIR/'annotation',
        scratch = 1000,
        mem = 4000,
        time = 235,
        qerrfile = OUTDIR/'annotation/{sample}/{sample}.prokka.qerr',
        qoutfile = OUTDIR/'annotation/{sample}/{sample}.prokka.qout'
    conda:
        'envs/annotate.yaml'
    threads:
        2
    shell:
        'prokka --outdir {params.outdir}/{params.locustag} '
        '--locustag {params.locustag} '
        '--compliant '
        '--prefix {params.locustag} {input.scaffolds} '
        '--force '


rule prokka_plasmid:
    input:
        scaffolds = OUTDIR/'plasmid/{sample}/scaffolds.fasta',
        marker = OUTDIR/'plasmid/{sample}/{sample}.spades.done'
    output:
        gff = OUTDIR/'plasmid/annotation/{sample}/{sample}.gff',
        gbk = OUTDIR/'plasmid/annotation/{sample}/{sample}.gbk',
        fna = OUTDIR/'plasmid/annotation/{sample}/{sample}.fna',
        faa = OUTDIR/'plasmid/annotation/{sample}/{sample}.faa',
        marker = touch(OUTDIR/'plasmid/annotation/{sample}/{sample}.prokka.done')
    params:
        locustag = '{sample}',
        outdir = OUTDIR/'plasmid/annotation',
        scratch = 1000,
        mem = 4000,
        time = 235,
        qerrfile = OUTDIR/'plasmid/annotation/{sample}/{sample}.prokka.qerr',
        qoutfile = OUTDIR/'plasmid/annotation/{sample}/{sample}.prokka.qout'
    conda:
        'envs/annotate.yaml'
    threads:
        2
    shell:
        'prokka --outdir {params.outdir}/{params.locustag} '
        '--locustag {params.locustag} '
        '--compliant '
        '--prefix {params.locustag} {input.scaffolds} '
        '--force '



# rule gunzip:
#     input: '{sample}.fasta.gz'
#     output: '{sample}.fasta'
#     params:
#         qerrfile = '{sample}.gzip.qerr',
#         qoutfile = '{sample}.gzip.qout',
#         scratch = 6000,
#         mem = 7700,
#         time = 1400
#
#     threads:
#         8
#     shell:
#         'gunzip {input}'


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

