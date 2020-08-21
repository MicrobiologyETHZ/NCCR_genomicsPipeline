from pathlib import Path
import sys

DATADIR = Path(config['dataDir'])
OUTDIR = Path(config['outDir'])
ADAPTERS = Path(config['adapters'])
PHIX = Path(config['phix'])
SFILE = Path(config['sampleFile'])

def get_merged(wildcards):
    if config['merged']== True:
        return '--pe-1-m ' + str(OUTDIR) + f'/merged_reads/{wildcards.subsample}/{wildcards.sample}.m.fq.gz '
    else:
        return  ''

READ_DIR = 'merged_reads' if config['merged'] else 'clean_reads'

# todo REMOVE DUPLICATED CODE
def get_subsamples(sample_file=SFILE):
    subsamples = []
    if Path(sample_file).is_file():
        subsamples = set(Path(sample_file).read_text().splitlines())
    if len(subsamples) == 0:
        exit(1)
    return subsamples



def get_samples(subsamples, suffix = '.fq.gz'):
    samples = []

    for s in subsamples:
        pattern = f'{s}*{suffix}'
        sample_path = list(Path(DATADIR).joinpath(s).rglob(pattern))
        if len(sample_path)<1:
            sys.exit(1)
        else:
            sample = str(sample_path[0].name).split('.')[0]
            samples.append(sample)
    return samples


SUBSAMPLES = get_subsamples()
SAMPLES = get_samples(SUBSAMPLES)



rule assemble:
        input: [OUTDIR/f'assembly/{sub}/{sam}.spades.done' for sub, sam in zip(SUBSAMPLES, SAMPLES)],
            [OUTDIR/f'assembly/plasmid/{sub}/{sam}.plasmid.spades.done' for sub, sam in zip(SUBSAMPLES, SAMPLES)]
# todo Incorporate merged_reads as an option

rule assemble_wga:
        input:
            fq1 = OUTDIR/READ_DIR/'{subsample}/{sample}.1.fq.gz',
            fq2 = OUTDIR/READ_DIR/'{subsample}/{sample}.2.fq.gz',
            s = OUTDIR/READ_DIR/'{subsample}/{sample}.s.fq.gz',
        output:
            marker = touch(OUTDIR/'assembly/{subsample}/{sample}.spades.done'),
        #marker_merge = touch('{path}/'+SPADES_FOLDER_NAME+'/{sample}.spades.merge.done')
        params:
            subsample = '{subsample}',
            outdir = OUTDIR/'assembly/',
            qerrfile = OUTDIR/'assembly/{subsample}/{sample}.spades.qerr',
            qoutfile = OUTDIR/'assembly/{subsample}/{sample}.spades.qout',
            merged = get_merged,
            scratch = 6000,
            mem = 7700,
            time = 1400
        log:
            log = OUTDIR/'assembly/{subsample}/{sample}.spades.log',
            command = OUTDIR/'assembly/{subsample}/{sample}.spades.command'
        conda:
            'envs/assemble.yaml'
        benchmark:
            OUTDIR/'assembly/{subsample}/{sample}.spades.benchmark'
        threads:
            8
        shell:
            "spades.py -t 4 --careful "
            " --pe1-1 {input.fq1} --pe1-2 {input.fq2} "
            "--pe1-s {input.s} {params.merged}"
            "-o {params.outdir}/{params.subsample} &> {log.log} "



rule assemble_plasmid:
        input:
            fq1 = OUTDIR/READ_DIR/'{subsample}/{sample}.1.fq.gz',
            fq2 = OUTDIR/READ_DIR/'{subsample}/{sample}.2.fq.gz',
            s = OUTDIR/READ_DIR/'{subsample}/{sample}.s.fq.gz',
        output:
            marker = touch(OUTDIR/'assembly/plasmid/{subsample}/{sample}.plasmid.spades.done'), #todo put these all in the OUTDIR
        params:
            subsample = '{subsample}',
            outdir = OUTDIR/'assembly/plasmid/',
            qerrfile = OUTDIR/'assembly/plasmid/{subsample}/{sample}.spades.qerr',
            qoutfile = OUTDIR/'assembly/plasmid/{subsample}/{sample}.spades.qout',
            merged = get_merged,
            scratch = 6000,
            mem = 7700,
            time = 1400
        log:
            log = OUTDIR/'assembly/plasmid/{subsample}/{sample}.spades.log',
            command = OUTDIR/'assembly/plasmid/{subsample}/{sample}.spades.command'
        conda:
            'envs/assemble.yaml'
        benchmark:
            OUTDIR/'assembly/plasmid/{subsample}/{sample}.spades.benchmark'
        threads:
            8
        shell:
            "spades.py -t 4 --careful "
            " --pe1-1 {input.fq1} --pe1-2 {input.fq2} "
            "--pe1-s {input.s} {params.merged} "
            "--plasmid "
            "-o {params.outdir}/{params.subsample} &> {log.log} "


##############################


# rule assembly_cleanup:
#     input:
#         marker = '{path}/'+SPADES_FOLDER_NAME+'/{sample}.spades.done',
#         marker_plasmid =''
#     output:
#         marker = touch('{path}/'+SPADES_FOLDER_NAME+'/{sample}.assembly_cleanup.done'),
#         scaffolds0 = '{path}/'+SPADES_FOLDER_NAME+'/{sample}.scaffolds.min0.fasta.gz',
#         scaffolds500 = '{path}/'+SPADES_FOLDER_NAME+'/{sample}.scaffolds.min500.fasta.gz',
#         scaffolds1000 = '{path}/'+SPADES_FOLDER_NAME+'/{sample}.scaffolds.min1000.fasta.gz',
#         contigs0 = '{path}/'+SPADES_FOLDER_NAME+'/{sample}.contigs.min0.fasta.gz',
#         contigs500 = '{path}/'+SPADES_FOLDER_NAME+'/{sample}.contigs.min500.fasta.gz',
#         contigs1000 = '{path}/'+SPADES_FOLDER_NAME+'/{sample}.contigs.min1000.fasta.gz',
#         stats = '{path}/'+SPADES_FOLDER_NAME+'/{sample}.assembly.stats',
#     params:
#         mem = 1000,
#         scratch = 1000,
#         time = 20,
#         qerrfile = '{path}/'+SPADES_FOLDER_NAME+'/{sample}.assembly_cleanup.qerr',
#         qoutfile = '{path}/'+SPADES_FOLDER_NAME+'/{sample}.assembly_cleanup.qout',
#         workfolder = '{path}/'+SPADES_FOLDER_NAME+'/',
#         sample = '{sample}'
#     log:
#         command = '{path}/'+SPADES_FOLDER_NAME+'/{sample}.assembly_cleanup.stats.command'
#     threads:
#         8
#     shell:
#         '''
#         #!/bin/bash
#
#         command="
#         if [[ -f "{params.workfolder}scaffolds.fasta" ]]; then
#             pigz -p {threads} {params.workfolder}scaffolds.fasta
#         fi
#         if [[ -f "{params.workfolder}contigs.fasta" ]]; then
#             pigz -p {threads} {params.workfolder}contigs.fasta
#         fi
#         if [[ -f "{params.workfolder}assembly_graph.fastg" ]]; then
#             pigz -p {threads} {params.workfolder}assembly_graph.fastg
#         fi
#         if [[ -f "{params.workfolder}assembly_graph_with_scaffolds.gfa" ]]; then
#             pigz -p {threads} {params.workfolder}assembly_graph_with_scaffolds.gfa
#         fi
#         if [[ -f "{params.workfolder}contigs.paths" ]]; then
#             pigz -p {threads} {params.workfolder}contigs.paths
#         fi
#         if [[ -f "{params.workfolder}scaffolds.paths" ]]; then
#             pigz -p {threads} {params.workfolder}scaffolds.paths
#         fi
#         if [[ -f "{params.workfolder}misc/broken_scaffolds.fasta" ]]; then
#             rm {params.workfolder}misc/broken_scaffolds.fasta
#         fi
#         if [[ -f "{params.workfolder}first_pe_contigs.fasta" ]]; then
#             rm {params.workfolder}first_pe_contigs.fasta
#         fi
#         if [[ -f "{params.workfolder}before_rr.fasta" ]]; then
#             rm {params.workfolder}before_rr.fasta
#         fi
#         python /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/PAN/GENERAL_METAGT_ANALYSIS_PAN/code/pipeline/contig_filter.py {params.sample} contigs {params.workfolder}contigs.fasta.gz {params.workfolder}
#         python /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/PAN/GENERAL_METAGT_ANALYSIS_PAN/code/pipeline/contig_filter.py {params.sample} scaffolds {params.workfolder}scaffolds.fasta.gz {params.workfolder}
#         pigz -p {threads} {params.workfolder}*min*fasta
#         pigz -p {threads} {params.workfolder}*hashes
#         assembly-stats -l 500 -t <(zcat {output.scaffolds500}) > {output.stats}
#         ";
#         echo "$command" > {log.command};
#         eval "$command"
#         '''
#
#
#
#

  # # Clean paired and unpaired files exists. Take all these files as input.
  #       cmdstring = spades.py + " " + ConfigSectionMap("spades", Config)['spades_parameters'] + " --pe1-1 " + forward_paired + " --pe1-2 " + reverse_paired + " --pe1-s " + forward_unpaired + " --pe1-s " + reverse_unpaired + " -o " + out_path + "spades_results"
  #       plasmid_cmdstring = ConfigSectionMap("spades", Config)['base_cmd'] + " " + ConfigSectionMap("spades", Config)['plasmid_spades_parameters'] + " --pe1-1 " + forward_paired + " --pe1-2 " + reverse_paired + " --pe1-s " + forward_unpaired + " --pe1-s " + reverse_unpaired + " -o " + out_path + "spades_plasmid_results"




