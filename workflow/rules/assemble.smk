import sys
from pathlib import Path


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
            scaffolds = OUTDIR/'assembly/{sample}/scaffolds.fasta',
            contigs = OUTDIR/'assembly/{sample}/contigs.fasta',
            min200 = OUTDIR/'assembly/{sample}/{sample}.scaffolds.min200.fasta'
        params:
            sample = lambda wildcards: f'{wildcards.sample}',
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
            'assembly'
        threads:
            16
        shell:
            "spades.py -t 4 --isolate "
            " --pe1-1 {input.fq1} --pe1-2 {input.fq2} "
            "--pe1-s {input.s} {params.merged}"
            "-o {params.outdir} &> {log.log};"
            "python ./scripts/contig_filter.py {params.sample} contigs {params.outdir}/contigs.fasta {params.outdir} ISO ; "
            "python ./scripts/contig_filter.py {params.sample} scaffolds {params.outdir}/scaffolds.fasta {params.outdir} ISO "


rule unicycler_short:
    input:
        fq1=OUTDIR / READ_DIR / '{sample}/{sample}.1.fq.gz',
        fq2=OUTDIR / READ_DIR / '{sample}/{sample}.2.fq.gz',
        s=OUTDIR / READ_DIR / '{sample}/{sample}.s.fq.gz',
    output:
        marker = touch(OUTDIR/f'unicycler/{config["unimode"]}'/'{sample}/{sample}.unicycler.done'),
        assembly = OUTDIR/f'unicycler/{config["unimode"]}'/'{sample}/assembly.fasta'
    params:
        outdir = lambda wildcards: OUTDIR/f'unicycler/{config["unimode"]}/{wildcards.sample}',
        qerrfile = lambda wildcards: OUTDIR/f'logs/unicycler/{wildcards.sample}.unicycler.{config["unimode"]}.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/unicycler/{wildcards.sample}.unicycler.{config["unimode"]}.qout',
        threads = 24,
        mode = config["unimode"],
        scratch = 6000,
        mem = 7700,
        time = 1400
    log:
        log = OUTDIR/'logs/unicycler/{sample}.unicycler.log',
    conda:
        'assembly'
    threads:
        32
    shell:
        "unicycler -1 {input.fq1} -2 {input.fq2} "
        "-s {input.s} --mode {params.mode} "
        "-t {params.threads} -o {params.outdir} &> {log.log} "


rule prokka:
    input: '{assembly}.fasta'
    output: marker = touch('{assembly}.prokka.done')
    params:
        locustag=lambda wildcards: Path(f'{wildcards.assembly}').parent.stem,
        outdir=lambda wildcards: Path(f'{wildcards.assembly}').parent,
        scratch=1000,
        mem=4000,
        time=235,
        qerrfile=lambda wildcards: str(OUTDIR / 'logs' / Path(f'{wildcards.assembly}').parent.stem) + '.prokka.qerr',
        qoutfile=lambda wildcards: str(OUTDIR / 'logs' / Path(f'{wildcards.assembly}').parent.stem) + '.prokka.qout'
    conda:
        'assembly'
    log:
        log='{assembly}.prokka.log',
    threads:
        8
    shell:
        'prokka --outdir {params.outdir} '
        '--locustag {params.locustag} '
        '--compliant '
        '--prefix {params.locustag} {input} '
        '--force &> {log.log} '


if config['assembler'] == 'spades': # todo rethink cleanup

    rule assembly_cleanup:
        input:
            marker = OUTDIR/'{assembly}/{sample}/{sample}.spades.done',
            scaffolds = OUTDIR/'{assembly}/{sample}/scaffolds.fasta',
            contigs = OUTDIR/'{assembly}/{sample}/contigs.fasta'
        output:
            marker = touch(OUTDIR/'{assembly}/{sample}/{sample}.assembly_cleanup.done'),
            scaffolds0 = OUTDIR/'{assembly}/{sample}/{sample}.scaffolds.min0.fasta.gz',
            scaffolds200 = OUTDIR/'{assembly}/{sample}/{sample}.scaffolds.min200.fasta.gz',
            scaffolds1000 = OUTDIR/'{assembly}/{sample}/{sample}.scaffolds.min1000.fasta.gz',
            contigs0 = OUTDIR/'{assembly}/{sample}/{sample}.contigs.min0.fasta.gz',
            contigs200 = OUTDIR/'{assembly}/{sample}/{sample}.contigs.min200.fasta.gz',
            contigs1000 = OUTDIR/'{assembly}/{sample}/{sample}.contigs.min1000.fasta.gz',
            stats = OUTDIR/'{assembly}/{sample}/{sample}.assembly.stats'
        params:
            mem = 1000,
            scratch = 1000,
            time = 20,
            qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.assembly}/{wildcards.sample}.assembly_cleanup.qerr',
            qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.assembly}/{wildcards.sample}.assembly_cleanup.qout',
            workfolder = lambda wildcards: OUTDIR/f'{wildcards.assembly}',
            sample = '{sample}',
        conda:
            'assembly'
        log:
            log = OUTDIR/'logs/{assembly}/{sample}.assembly_cleanup.stats.log'
        threads:
            8
        shell:
            '''
            #!/bin/bash
    
            command="
            if [[ -f "{params.workfolder}/{params.sample}/scaffolds.fasta" ]]; then
                pigz   -p {threads} {params.workfolder}/{params.sample}/scaffolds.fasta
            fi
            if [[ -f "{params.workfolder}/{params.sample}/contigs.fasta" ]]; then
                pigz  -p {threads} {params.workfolder}/{params.sample}/contigs.fasta
            fi
            if [[ -f "{params.workfolder}/{params.sample}/assembly_graph.fastg" ]]; then
                pigz  -p {threads} {params.workfolder}/{params.sample}/assembly_graph.fastg
            fi
            if [[ -f "{params.workfolder}/{params.sample}/assembly_graph_with_scaffolds.gfa" ]]; then
                pigz  -p {threads} {params.workfolder}/{params.sample}/assembly_graph_with_scaffolds.gfa
            fi
            if [[ -f "{params.workfolder}/{params.sample}/contigs.paths" ]]; then
                pigz  -p {threads} {params.workfolder}/{params.sample}/contigs.paths
            fi
            if [[ -f "{params.workfolder}/{params.sample}/scaffolds.paths" ]]; then
                pigz  -p {threads} {params.workfolder}/{params.sample}/scaffolds.paths
            fi
            if [[ -f "{params.workfolder}/{params.sample}/misc/broken_scaffolds.fasta" ]]; then
                rm {params.workfolder}/{params.sample}/misc/broken_scaffolds.fasta
            fi
            if [[ -f "{params.workfolder}/{params.sample}/first_pe_contigs.fasta" ]]; then
                rm {params.workfolder}/{params.sample}/first_pe_contigs.fasta
            fi
            if [[ -f "{params.workfolder}/{params.sample}/before_rr.fasta" ]]; then
                rm {params.workfolder}/{params.sample}/before_rr.fasta
            fi
            python ./scripts/contig_filter.py {params.sample} contigs {params.workfolder}/{params.sample}/contigs.fasta.gz {params.workfolder}/{params.sample}
            python ./scripts/contig_filter.py {params.sample} scaffolds {params.workfolder}/{params.sample}/scaffolds.fasta.gz {params.workfolder}/{params.sample}
            pigz -f -p {threads} {params.workfolder}/{params.sample}/*min*fasta
            pigz -f -p {threads} {params.workfolder}/{params.sample}/*hashes
            assembly-stats -l 200 -t <(zcat {params.workfolder}/{params.sample}/{params.sample}.scaffolds.min200.fasta.gz) > {params.workfolder}/{params.sample}/{params.sample}.assembly.stats
            &> {log.log}
            ";
            eval "$command"
            '''

elif config['assembler'] == 'unicycler':
    pass

else:
    sys.exit(1)







