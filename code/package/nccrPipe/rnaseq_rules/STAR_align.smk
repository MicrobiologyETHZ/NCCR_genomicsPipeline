from pathlib import Path

rule STAR_index:
    input: config['refGenome']
    output: touch(OUTDIR/f'{config["ProjectName"]}.index.done')
    params:
        qerrfile = OUTDIR/'logs/STAR.index.qerr',
        qoutfile = OUTDIR/'logs/STAR.index.qout',
        genomeDir = config["genomeDir"],
        annotation = config["refAnn"],
        overhang = config["overhang"],
        threads = 32,
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/STAR_featurecounts.yaml'
    log: OUTDIR/'logs/STAR.index.log'
    threads:
        32
    shell: "STAR --runThreadN {params.threads} "
           "--runMode genomeGenerate "
           "--genomeFastaFiles {input} "
           "--genomeDir {params.genomeDir} "
           "--sjdbGTFfile {params.annotation} "
           "--sjdbOverhang {params.overhang} "
           #"--sjdbGTFtagExonParentTranscript Parent " # For GFF3 annotations


rule STAR_align:
    input: index_done = OUTDIR/f'{config["ProjectName"]}.index.done',
        fq1 = OUTDIR/'clean_reads/{sample}/{sample}.1.fq.gz',
        fq2= OUTDIR/'clean_reads/{sample}/{sample}.2.fq.gz'
    output: marker = touch(OUTDIR/'bam/{sample}/{sample}.done'),
         bam = OUTDIR/'bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam'
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.star_align.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.star_align.qout',
        genomeDir = config["genomeDir"],
        annotation = config["refAnn"],
        overhang = config["overhang"],
        maxIntron=config['maxIntron'],
        prefix = lambda wildcards: OUTDIR/f'bam/{wildcards.sample}/{wildcards.sample}_',
        threads = 8,
        scratch = 6000,
        mem = 8000,
        time = 1400
    conda:
        'envs/STAR_featurecounts.yaml'
    log: OUTDIR/'logs/{sample}.star_align.log'
    threads:
        32
    shell: "STAR --runThreadN {params.threads} "
           "--readFilesIn {input.fq1} {input.fq2} "
           "--readFilesCommand gunzip -c " #format for paired end
           "--genomeDir {params.genomeDir} "
           "--sjdbGTFfile {params.annotation} "
           "--sjdbOverhang {params.overhang} "
           "--outFileNamePrefix {params.prefix} "
           "--outSAMtype BAM SortedByCoordinate "
           "--outSAMunmapped Within "
           "--outSAMattributes Standard "
           "--alignIntronMax {params.maxIntron} "
           "--quantMode GeneCounts &> {log} "


rule featureCounts:
    input: bam = OUTDIR/'bam/{sample}/{sample}_Aligned.sortedByCoord.out.bam',
    output: OUTDIR/'counts/{sample}/{sample}.count.txt'
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.featurecounts.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}.featurecounts.qout',
        annotation = config["refAnn"],
        attribute = config["attribute"],
        strand = config["strand"],
        threads = 32,
        scratch = 6000,
        mem = 8000,
        time = 1400
    conda:
        'envs/STAR_featurecounts.yaml'
    log: OUTDIR/'logs/{sample}.featurecounts.log'
    threads:
        32
    shell:
        "featureCounts -p -T {params.threads} " 
        "-a {params.annotation} -o {output} " 
        "-g {params.attribute} {input.bam} -s {params.strand} &> {log}"

#TETSETE
# testtest test

