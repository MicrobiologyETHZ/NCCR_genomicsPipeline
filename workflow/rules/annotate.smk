#!/bin/bash

# rule zip:
#     input: [OUTDIR/f'assembly/{sample}/scaffolds.fasta.gz' for sample in SUBSAMPLES]
#
from pathlib import Path
import pandas as pd

"""
PGAP — NCBI's Prokaryotic Genome Annotation Pipeline, run over an arbitrary set
of genomes rather than samples.csv. Genomes come from a samplesheet generated
by `nccrPipe pgap-samples` (see workflow/scripts/pgap_samplesheet.py), with
columns name,fasta,taxon — taxon is per-genome because PGAP's `-s` can't be
inferred from a filename.

Config:
    pgap_samples: path/to/pgap_samples.csv   # optional; PGAP_SAMPLES is {} if unset
    pgap:
      pgap_dir: /path/to/pgap_install         # required if pgap_samples is set
      cache: /path/to/pgap_install/cache      # optional, defaults to pgap_dir/cache
"""
_pgap_samples_file = config.get('pgap_samples', '')
if _pgap_samples_file:
    _pgap_df = pd.read_csv(_resolve(_pgap_samples_file), comment='#')
    PGAP_SAMPLES = {
        row['name']: {'fasta': str(_resolve(row['fasta'])), 'taxon': row['taxon']}
        for _, row in _pgap_df.iterrows()
    }
else:
    PGAP_SAMPLES = {}

_pgap_cfg = config.get('pgap', {})
PGAP_DIR = str(_resolve(_pgap_cfg['pgap_dir'])) if _pgap_cfg.get('pgap_dir') else ''
PGAP_CACHE = (str(_resolve(_pgap_cfg['cache'])) if _pgap_cfg.get('cache')
              else f'{PGAP_DIR}/cache')

if PGAP_SAMPLES and not PGAP_DIR:
    raise ValueError(
        "pgap_samples is set but pgap.pgap_dir is missing from the config.\n"
        "Add:\n"
        "  pgap:\n"
        "    pgap_dir: /path/to/pgap_install\n"
        "    # cache: /path/to/pgap_install/cache   # optional, defaults to pgap_dir/cache"
    )

wildcard_constraints:
    name = r'[A-Za-z0-9._\-]+'


rule pgap:
    """Annotate one genome with PGAP, running inside its own Apptainer container.

    Four things that matter here, all learned the hard way:
      1. `unset SLURM_CPUS_PER_TASK NSLOTS` — otherwise PGAP passes --cpus to
         Apptainer and hits a cgroup-v2 crash. threads: 8 above still reserves
         8 cores via the SLURM submission (params.mem is per-cpu); this only
         hides the count from PGAP itself.
      2. No `conda:` directive — PGAP needs a clean host interpreter; every
         tool it runs lives inside its own container.
      3. `directory()` output, not a file path inside it — Snakemake
         pre-creates the parent dirs of file outputs, and PGAP refuses a
         pre-existing -o dir. `rm -rf` first clears leftovers from a failed run.
      4. Apptainer must already be on PATH (confirmed true on this cluster).
    """
    input:
        fasta = lambda wc: PGAP_SAMPLES[wc.name]['fasta']
    output:
        outdir = directory(OUTDIR/'pgap/{name}')
    params:
        taxon = lambda wc: PGAP_SAMPLES[wc.name]['taxon'],
        pgap = f'{PGAP_DIR}/pgap.py',
        cache = PGAP_CACHE,
        scratch = 20000,
        mem = 32000,
        time = 300,
        qerrfile = lambda wc: OUTDIR/f'logs/pgap/{wc.name}.qerr',
        qoutfile = lambda wc: OUTDIR/f'logs/pgap/{wc.name}.qout'
    log:
        log = OUTDIR/'logs/pgap/{name}.log'
    threads:
        8
    shell:
        r"""
        rm -rf {output.outdir}
        unset SLURM_CPUS_PER_TASK NSLOTS

        export SSL_CERT_FILE=/etc/pki/tls/certs/ca-bundle.crt
        export REQUESTS_CA_BUNDLE=/etc/pki/tls/certs/ca-bundle.crt
        export PGAP_INPUT_DIR={params.cache}

        {{ echo "PGAP cache version: $(cat {params.cache}/VERSION 2>/dev/null)"; {params.pgap} --version; }} > {log.log} 2>&1

        {params.pgap} -n -g {input.fasta} -s "{params.taxon}" \
            -D apptainer --auto-correct-tax -o {output.outdir} &>> {log.log}
        """

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
db = config.get('database', 'eggnog')
if db == 'eggnog':
    rule emapper:
        input: faa = OUTDIR/"{assembly}/{sample}/{sample}.faa"
        #input: faa = OUTDIR/"{assembly}/{sample}/prokka/{sample}.faa"
        output: marker = touch(OUTDIR/"{assembly}/{sample}/eggnog/{sample}.eggnog.done")
        params:
            sample = "{sample}",
            outdir = lambda wildcards: OUTDIR/f'{wildcards.assembly}/{wildcards.sample}/eggnog',
            dataDir = config.get('eggnog_db', ''),
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

elif db == 'kegg':
    rule kegg:
        input: faa = OUTDIR/'{assembly}/{sample}/{sample}.faa'
        output: touch(OUTDIR/'{assembly}/{sample}/kegg/{sample}.kegg.done')
        params:
            outdir = lambda wildcards: OUTDIR/f'{wildcards.assembly}/{wildcards.sample}/kegg/', # needs to exist before starting the rule
            prefix = lambda wildcards: OUTDIR/f'{wildcards.assembly}/{wildcards.sample}/kegg/{wildcards.sample}',
            dataDir = f"/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/PAN/GENOMES_COLLECTION_PAN/data/resources/soft/kegg_annotation/apr2022/kegg_db/kegg/kegg_for_prokka-with_ko.dmnd", # Path to the (processed) kegg database
            koDir = f"/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/PAN/GENOMES_COLLECTION_PAN/data/resources/soft/kegg_annotation/apr2022/kegg_db/kegg/ko", # Path to the ko directory containing the mapping files
            scratch = 1000,
            mem = 4000,
            time = 235,
            qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.assembly}/{wildcards.sample}/kegg/{wildcards.sample}.kegg.qerr',
            qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.assembly}/{wildcards.sample}/kegg/{wildcards.sample}.kegg.qout'
        conda:
            'KEGG'
        log:
            log =  OUTDIR/'logs/{assembly}/{sample}/kegg/{sample}.kegg.log'
        threads:
            16
        shell:
            'python /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/PAN/GENOMES_COLLECTION_PAN/data/resources/soft/kegg_annotation/annotate_kegg_v2.py '
            '-r {params.outdir} -q {input} -o {params.prefix} -t {threads} '
            '-d {params.dataDir} '
            '-k {params.koDir}'

else:
    print("Error: Please choose a valid database for functional annotation. The options are eggnog or kegg.")
    sys.exit(1)


rule gapseq_find:
    input:
        faa = OUTDIR/"gapseq_data/{genome_prefix}.faa"
    output:
        #pathways = OUTDIR/"gapseq_output/{genome_prefix}/{genome_prefix}-all-Pathways.tbl",
        tcs = OUTDIR/"gapseq_output/{genome_prefix}/{genome_prefix}-Transporter.tbl",
        #reactions = "gapseq_output/{genome}/{genome}-all-Reactions.tbl",
        #transporter = "gapseq_output/{genome}/{genome}-Transporter.tbl"
    params:
        outdir = lambda wildcards: OUTDIR/f"gapseq_output/{wildcards.genome_prefix}",
        genome_id = lambda wildcards: wildcards.genome_prefix,
        scratch = 1000,
        mem = 4000,
        time = 235,
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.genome_prefix}.gapseq.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.genome_prefix}.gapseq.qout'
    conda:
        "gapseq"  
    log:
        OUTDIR/"logs/gapseq_find/{genome_prefix}.log"
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        
        #gapseq find -p all -b 200 -t Bacteria -m Bacteria \
        #   {input.faa} &> {log}
        gapseq find-transport {input.faa}  &> {log}

        # Rename outputs to include genome name
        #mv *-all-Pathways.tbl {params.genome_id}-all-Pathways.tbl || true
        #mv *-all-Reactions.tbl {params.genome_id}-all-Reactions.tbl || true
        mv *Transporter.tbl {params.genome_id}-Transporter.tbl || true
        """



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

