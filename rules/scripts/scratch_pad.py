from pathlib import Path

def get_samples(sample_file=Path('/nfs/home/ansintsova/esbl_strains/samples.txt'),
                dataDir = Path('/nfs/home/ansintsova/esbl_strains/raw')):
    samples = []
    raw_seqfiles = []
    if sample_file.is_file():
        samples = set(sample_file.read_text().splitlines())
    if len(samples) == 0:
        exit(1)
    for sample in samples:
        print(dataDir.joinpath(sample))
        for raw_file in dataDir.joinpath(sample).rglob('*.fq.gz'):

            raw_seqfiles.append(raw_file)
    return samples, raw_seqfiles

print(get_samples())

rule assembly_cleanup:
    input:
        marker = OUTDIR/'assembly/{subsample}/{sample}.spades.done',
        marker_plasmid =OUTDIR/'assembly/plasmid/{subsample}/{sample}.plasmid.spades.done'
    output:
        marker = touch(OUTDIR/'assembly/{subsample}/{sample}.assembly_cleanup.done'),
        scaffolds0 = OUTDIR/'assembly/{subsample}/{sample}.scaffolds.min0.fasta.gz',
        scaffolds500 = OUTDIR/'assembly/{subsample}/{sample}.scaffolds.min500.fasta.gz',
        scaffolds1000 = OUTDIR/'assembly/{subsample}/{sample}.scaffolds.min1000.fasta.gz',
        contigs0 = OUTDIR/'assembly/{subsample}/{sample}.contigs.min0.fasta.gz',
        contigs500 = OUTDIR/'assembly/{subsample}/{sample}.contigs.min500.fasta.gz',
        contigs1000 = OUTDIR/'assembly/{subsample}/{sample}.contigs.min1000.fasta.gz',
        stats = OUTDIR/'assembly/{subsample}/{sample}.assembly.stats',
    params:
        mem = 1000,
        scratch = 1000,
        time = 20,
        qerrfile = OUTDIR/'assembly/{subsample}/{sample}.assembly_cleanup.qerr',
        qoutfile = OUTDIR/'assembly/{subsample}/{sample}.assembly_cleanup.qout',
        workfolder = OUTDIR/'assembly/{subsample}/',
        sample = '{sample}'
    log:
        command = OUTDIR/'assembly/{subsample}/{sample}.assembly_cleanup.stats.command'
    threads:
        8
    shell:
        '''
        #!/bin/bash

        command="
        if [[ -f "{params.workfolder}scaffolds.fasta" ]]; then
            pigz -p {threads} {params.workfolder}scaffolds.fasta
        fi
        if [[ -f "{params.workfolder}contigs.fasta" ]]; then
            pigz -p {threads} {params.workfolder}contigs.fasta
        fi
        if [[ -f "{params.workfolder}assembly_graph.fastg" ]]; then
            pigz -p {threads} {params.workfolder}assembly_graph.fastg
        fi
        if [[ -f "{params.workfolder}assembly_graph_with_scaffolds.gfa" ]]; then
            pigz -p {threads} {params.workfolder}assembly_graph_with_scaffolds.gfa
        fi
        if [[ -f "{params.workfolder}contigs.paths" ]]; then
            pigz -p {threads} {params.workfolder}contigs.paths
        fi
        if [[ -f "{params.workfolder}scaffolds.paths" ]]; then
            pigz -p {threads} {params.workfolder}scaffolds.paths
        fi
        if [[ -f "{params.workfolder}misc/broken_scaffolds.fasta" ]]; then
            rm {params.workfolder}misc/broken_scaffolds.fasta
        fi
        if [[ -f "{params.workfolder}first_pe_contigs.fasta" ]]; then
            rm {params.workfolder}first_pe_contigs.fasta
        fi
        if [[ -f "{params.workfolder}before_rr.fasta" ]]; then
            rm {params.workfolder}before_rr.fasta
        fi
        python scripts/contig_filter.py {params.sample} contigs {params.workfolder}contigs.fasta.gz {params.workfolder}
        python scripts/contig_filter.py {params.sample} scaffolds {params.workfolder}scaffolds.fasta.gz {params.workfolder}
        pigz -p {threads} {params.workfolder}*min*fasta
        pigz -p {threads} {params.workfolder}*hashes
        assembly-stats -l 500 -t <(zcat {output.scaffolds500}) > {output.stats}
        ";
        echo "$command" > {log.command};
        eval "$command"
        '''