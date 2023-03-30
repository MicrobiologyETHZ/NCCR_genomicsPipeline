#
# rule zip:
#     input: [OUTDIR/f'assembly/{sample}/scaffolds.fasta.gz' for sample in SUBSAMPLES]
#
from pathlib import Path

#
# rule gunzipAnn:
#     input:
#         OUTDIR/'{assembly}/{sample1}/{sample1}.scaffolds.min0.fasta.gz'
#     output:
#         OUTDIR/'{assembly}/{sample1}/{sample1}.scaffolds.min0.fasta'
#     params:
#         qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.assembly}/{wildcards.sample1}.gzip.qerr',
#         qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.assembly}/{wildcards.sample1}.gzip.qout',
#         scratch = 6000,
#         mem = 7700,
#         time = 1400
#     conda:
#         'envs/annotate.yaml'
#     threads:
#         8
#     shell:
#         "gunzip {input} "


rule emapper:
    input: faa = OUTDIR/"{assembly}/{sample}/{sample}.faa"
    #input: faa = OUTDIR/"{assembly}/{sample}/prokka/{sample}.faa"
    output: marker = touch(OUTDIR/"{assembly}/{sample}/eggnog/{sample}.eggnog.done")
    params:
        sample = "{sample}",
        outdir = lambda wildcards: OUTDIR/f'{wildcards.assembly}/{wildcards.sample}/eggnog',
        dataDir = "/science/ansintsova/eggnog-data/", # todo put this into config
        scratch = 1000,
        mem = 4000,
        time = 235,
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.assembly}/{wildcards.sample}/eggnog/{wildcards.sample}.emapper.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.assembly}/{wildcards.sample}/eggnog/{wildcards.sample}.emapper.qout'
    conda:
        'emapper'
    log:
        log = OUTDIR/'logs/{assembly}/{sample}/eggnog/{sample}.emapper.log'
    threads:
        16
    shell:
        'emapper.py -i {input.faa} --output_dir {params.outdir} --output {params.sample} '
        '--cpu 16 --temp_dir {params.outdir} '
        ' -m diamond --data_dir {params.dataDir} &> {log.log} '


# todo add KEGG annotation


#
# rule prokka_plasmid:
#     input:
#         scaffolds = OUTDIR/'plasmid/{sample}/scaffolds.fasta',
#         marker = OUTDIR/'plasmid/{sample}/{sample}.spades.done'
#     output:
#         gff = OUTDIR/'plasmid/annotation/{sample}/{sample}.gff',
#         gbk = OUTDIR/'plasmid/annotation/{sample}/{sample}.gbk',
#         fna = OUTDIR/'plasmid/annotation/{sample}/{sample}.fna',
#         faa = OUTDIR/'plasmid/annotation/{sample}/{sample}.faa',
#         marker = touch(OUTDIR/'plasmid/annotation/{sample}/{sample}.prokka.done')
#     params:
#         locustag = '{sample}',
#         outdir = OUTDIR/'plasmid/annotation',
#         scratch = 1000,
#         mem = 4000,
#         time = 235,
#         qerrfile = OUTDIR/'plasmid/annotation/{sample}/{sample}.prokka.qerr',
#         qoutfile = OUTDIR/'plasmid/annotation/{sample}/{sample}.prokka.qout'
#
#     conda:
#         'envs/annotate.yaml'
#     threads:
#         2
#     shell:
#         'prokka --outdir {params.outdir}/{params.locustag} '
#         '--locustag {params.locustag} '
#         '--compliant '
#         '--prefix {params.locustag} {input.scaffolds} '
#         '--force '




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

