

def get_merged(wildcards):
    if config['merged']== True:
        return '--pe1-m ' + str(OUTDIR) + f'/merged_reads/{wildcards.sample}/{wildcards.sample}.m.fq.gz '
    else:
        return  ''

READ_DIR = 'merged_reads' if config['merged'] else 'clean_reads'

rule assemble_wga:
        input:
            fq1 = OUTDIR/READ_DIR/'{sample}/{sample}.1.fq.gz',
            fq2 = OUTDIR/READ_DIR/'{sample}/{sample}.2.fq.gz',
            s = OUTDIR/READ_DIR/'{sample}/{sample}.s.fq.gz',
        output:
            marker = touch(OUTDIR/'assembly/{sample}/{sample}.spades.done'),
            #scaffolds = OUTDIR/'assembly/{sample}/scaffolds.fasta',
            #contigs = OUTDIR/'assembly/{sample}/contigs.fasta'
        params:
            outdir = lambda wildcards: OUTDIR/f'assembly/{wildcards.sample}',
            qerrfile = lambda wildcards: OUTDIR/f'logs/assembly/{wildcards.sample}.spades.qerr',
            qoutfile = lambda wildcards: OUTDIR/f'logs/assembly/{wildcards.sample}.spades.qout',
            merged = get_merged,
            scratch = 6000,
            mem = 7700,
            time = 1400
        log:
            log = OUTDIR/'logs/assembly/{sample}.spades.log',
        conda:
            'envs/assemble.yaml'
        benchmark:
            OUTDIR/'assembly/{sample}/{sample}.spades.benchmark'
        threads:
            8
        shell:
            "spades.py -t 4 --isolate "
            " --pe1-1 {input.fq1} --pe1-2 {input.fq2} "
            "--pe1-s {input.s} {params.merged}"
            "-o {params.outdir} &> {log.log} "



rule assemble_plasmid:
        input:
            fq1 = OUTDIR/READ_DIR/'{sample}/{sample}.1.fq.gz',
            fq2 = OUTDIR/READ_DIR/'{sample}/{sample}.2.fq.gz',
            s = OUTDIR/READ_DIR/'{sample}/{sample}.s.fq.gz',
        output:
            marker = touch(OUTDIR/'plasmid/{sample}/{sample}.spades.done'), #todo put these all in the OUTDIR
            #scaffolds = OUTDIR/'plasmid/{sample}/scaffolds.fasta',
            #contigs = OUTDIR/'plasmid/{sample}/contigs.fasta'
        params:
            outdir = lambda wildcards: OUTDIR/f'assembly/{wildcards.sample}',
            qerrfile = lambda wildcards: OUTDIR/f'logs/plasmid/{wildcards.sample}.plasmid.spades.qerr',
            qoutfile = lambda wildcards: OUTDIR/f'logs/plasmid/{wildcards.sample}.plasmid.spades.qout',
            merged = get_merged,
            scratch = 6000,
            mem = 7700,
            time = 1400
        log:
            log = OUTDIR/'logs/plasmid/{sample}.plasmid.spades.log',
        conda:
            'envs/assemble.yaml'
        benchmark:
            OUTDIR/'plasmid/{sample}/{sample}.spades.benchmark'
        threads:
            8
        shell:
            "spades.py -t 4 --careful --plasmid"
            " --pe1-1 {input.fq1} --pe1-2 {input.fq2} "
            "--pe1-s {input.s} {params.merged} "
            "-o {params.outdir} &> {log.log} "


##############################
rule clean:
    input:
        [OUTDIR/f'assembly/{sub}/{sub}.assembly_cleanup.done' for sub in SUBSAMPLES],
        [OUTDIR/f'plasmid/{sub}/{sub}.assembly_cleanup.done' for sub in SUBSAMPLES]



rule assembly_cleanup:
    input:
        marker = OUTDIR/'{assembly}/{sample}.spades.done',
        scaffolds = OUTDIR/'{assembly}/scaffolds.fasta',
        contigs = OUTDIR/'{assembly}/contigs.fasta'
    output:
        marker = touch(OUTDIR/'{assembly}/{sample}.assembly_cleanup.done'),
        scaffolds0 = OUTDIR/'{assembly}/{sample}.scaffolds.min0.fasta.gz',
        scaffolds500 = OUTDIR/'{assembly}/{sample}.scaffolds.min500.fasta.gz',
        scaffolds1000 = OUTDIR/'{assembly}/{sample}.scaffolds.min1000.fasta.gz',
        contigs0 = OUTDIR/'{assembly}/{sample}.contigs.min0.fasta.gz',
        contigs500 = OUTDIR/'{assembly}/{sample}.contigs.min500.fasta.gz',
        contigs1000 = OUTDIR/'{assembly}/{sample}.contigs.min1000.fasta.gz',
        stats = OUTDIR/'{assembly}/{sample}.assembly.stats'
    params:
        mem = 1000,
        scratch = 1000,
        time = 20,
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.assembly}/{wildcards.sample}.assembly_cleanup.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.assembly}/{wildcards.sample}.assembly_cleanup.qout',
        workfolder = lambda wildcards: OUTDIR/f'{wildcards.assembly}',
        sample = '{sample}',
    conda:
        'envs/assembly_cleanup.yaml'
    log:
        log = OUTDIR/'logs/{assembly}/{sample}.assembly_cleanup.stats.log',
        command = OUTDIR/'logs/{assembly}/{sample}.assembly_cleanup.stats.command'
    threads:
        8
    shell:
        '''
        #!/bin/bash

        command="
        if [[ -f "{params.workfolder}/scaffolds.fasta" ]]; then
            pigz   -p {threads} {params.workfolder}/scaffolds.fasta
        fi
        if [[ -f "{params.workfolder}/contigs.fasta" ]]; then
            pigz  -p {threads} {params.workfolder}/contigs.fasta
        fi
        if [[ -f "{params.workfolder}/assembly_graph.fastg" ]]; then
            pigz  -p {threads} {params.workfolder}/assembly_graph.fastg
        fi
        if [[ -f "{params.workfolder}/assembly_graph_with_scaffolds.gfa" ]]; then
            pigz  -p {threads} {params.workfolder}/assembly_graph_with_scaffolds.gfa
        fi
        if [[ -f "{params.workfolder}/contigs.paths" ]]; then
            pigz  -p {threads} {params.workfolder}/contigs.paths
        fi
        if [[ -f "{params.workfolder}/scaffolds.paths" ]]; then
            pigz  -p {threads} {params.workfolder}/scaffolds.paths
        fi
        if [[ -f "{params.workfolder}/misc/broken_scaffolds.fasta" ]]; then
            rm {params.workfolder}/misc/broken_scaffolds.fasta
        fi
        if [[ -f "{params.workfolder}/first_pe_contigs.fasta" ]]; then
            rm {params.workfolder}/first_pe_contigs.fasta
        fi
        if [[ -f "{params.workfolder}/before_rr.fasta" ]]; then
            rm {params.workfolder}/before_rr.fasta
        fi
        python /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/PAN/GENERAL_METAGT_ANALYSIS_PAN/code/pipeline/contig_filter.py {params.sample} contigs {params.workfolder}/contigs.fasta.gz {params.workfolder}
        python /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/PAN/GENERAL_METAGT_ANALYSIS_PAN/code/pipeline/contig_filter.py {params.sample} scaffolds {params.workfolder}/scaffolds.fasta.gz {params.workfolder}
        pigz -f -p {threads} {params.workfolder}/*min*fasta
        pigz -f -p {threads} {params.workfolder}/*hashes
        assembly-stats -l 500 -t <(zcat {params.workfolder}/{params.sample}.scaffolds.min500.fasta.gz) > {params.workfolder}/{params.sample}.assembly.stats
        &> {log.log}
        ";
        echo "$command" > {log.command};
        eval "$command"
#         '''










