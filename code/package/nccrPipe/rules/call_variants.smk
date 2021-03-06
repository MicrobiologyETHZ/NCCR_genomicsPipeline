# if config['reference'].endswith('.fasta') or config['reference'].endswith('.fasta.gz') or config['reference'].endswith('.fna.gz'):
#     REFPATH = Path(config['reference'])
# elif config['scaffold'] == 'min0':
#     REFPATH = OUTDIR/f'assembly/{config["reference"]}/{config["reference"]}.scaffolds.min0.fasta'
#
#     # rule gunzip:
#     #     input: str(REFPATH) + ".gz"
#     #     output: REFPATH
#     #     shell: 'gunzip {input}'
#
# else:
#     REFPATH = OUTDIR/f'assembly/{config["reference"]}/scaffolds.fasta'


rule remove_duplicates:
    input: bam = OUTDIR/'{bam}/{sample}/{sample}.bam',
        inmar = OUTDIR/'{bam}/{sample}/{sample}.refbwa.done'
    output: bam = OUTDIR/'{bam}/{sample}/{sample}.rmdup.bam',
        marker = touch(OUTDIR/'{bam}/{sample}/{sample}.rmdup.done')
    params:
        qerrfile = lambda wildcards: OUTDIR/f'{wildcards.bam}/{wildcards.sample}/{wildcards.sample}.markdup.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'{wildcards.bam}/{wildcards.sample}/{wildcards.sample}.markdup.qout',
        scratch = 6000,
        mem = 10000,
        time = 1400,
        ram = config['ram'],
        tmpdir = config['tmpdir'],
        metrics = lambda wildcards: OUTDIR/f'{wildcards.bam}/{wildcards.sample}/{wildcards.sample}.picard.metrics'
    conda:
        'envs/call_variants.yaml'
    log:
        log = OUTDIR/'logs/{bam}/{sample}.removeDuplicates.log'
    threads:
        8
    shell: "gatk --java-options '-Xmx{params.ram}G' "
           "MarkDuplicates --TMP_DIR {params.tmpdir} "
           "--CREATE_INDEX true --REMOVE_DUPLICATES true "
           "-O {output.bam} -I {input.bam} -M {params.metrics} &> {log.log} "

# recalibrate quality?

rule bcf_call:
    input:
        scaf = REFPATH,
        bam = OUTDIR/"ref_bams/{sample}/{sample}.rmdup.bam"
    output:
        vcf = OUTDIR/'VCF/{sample}/{sample}.vcf',
        marker = touch(OUTDIR/'VCF/{sample}/{sample}.mpileup.done')
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}.mpileup.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}.mpileup.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    log:
        log = OUTDIR/'logs/VCF/{sample}.mpileup.log'
    conda:
        'envs/call_variants.yaml'
    threads:
        32
    shell:
        "bcftools mpileup -Ou -q30 -d3000 -f {input.scaf}  {input.bam} | " 
        "bcftools call -mv -o {output.vcf} 2> {log.log}"


rule sampile:
    input:
        scaf = REFPATH,
        bam = OUTDIR/"ref_bams/{sample}/{sample}.rmdup.bam"
    output:
        mpile = OUTDIR/'VCF/{sample}/{sample}.mpileup',
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}.smpileup.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}.smpileup.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    log:
        log = OUTDIR/'logs/VCF/{sample}.smpileup.log'
    conda:
        'envs/align.yaml'
    threads:
        32
    shell:
        "samtools mpileup -O -s -f {input.scaf}  {input.bam} -o {output.mpile} --output-extra QNAME"

# #
rule bcf_filter:
    input: vcf = OUTDIR/'VCF/{sample}/{sample}.vcf',
        m = OUTDIR/'VCF/{sample}/{sample}.mpileup.done'
    output: fvcf = OUTDIR/'VCF/{sample}/{sample}.filtered.vcf',
        marker = touch(OUTDIR/'VCF/{sample}/{sample}.vcf.done')
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}.bcf.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}.bcf.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    log:
        log = OUTDIR/'logs/VCF/{sample}.bcf.log'
    conda:
        'envs/call_variants.yaml'
    threads:
        8
    shell:
        "bcftools filter -Ov -sLowQual -g5 -G10 -e 'QUAL<100 || DP4[2]<10 || DP4[3]<10 ||  MQ<60' {input.vcf} | "
        "bcftools query  -i'FILTER=\"PASS\"' -f '%LINE' -o {output.fvcf} 2> {log.log}"


rule bcf_filter_isolate:
    input: OUTDIR/'VCF/{sample}/{sample}.vcf',

    output: fvcf = OUTDIR/'VCF/{sample}/{sample}.isolate.filtered.vcf',
        marker = touch(OUTDIR/'VCF/{sample}/{sample}.isolate.vcf.done')
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}.isolate.bcf.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}.isolate.bcf.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/call_variants.yaml'
    log:
        log = OUTDIR/'logs/VCF/{sample}.isolate.bcf.log'
    threads:
        8
    shell:
        "bcftools filter -Ov -sLowQual -g5 -G10 -e 'QUAL<10 ||  DP4[2]<10 || DP4[3]<10 ||(DP4[2] + DP4[3])/sum(DP4) < 0.9 ||  MQ<50' {input} | "
        "bcftools query  -i'FILTER=\"PASS\"' -f '%LINE' -o {output.fvcf} 2> {log.log}"



rule bcf_filter3:
    input: OUTDIR/'VCF/{sample}/{sample}.vcf',

    output: fvcf = OUTDIR/'VCF/{sample}/{sample}.f3.filtered.vcf',
        marker = touch(OUTDIR/'VCF/{sample}/{sample}.f3.vcf.done')
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}.f3.bcf.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}.f3.bcf.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/call_variants.yaml'
    log:
        log = OUTDIR/'logs/VCF/{sample}.f3.bcf.log'
    threads:
        8
    shell:
        "bcftools filter -Ov -sLowQual -g5 -G10 -e 'QUAL<10 ||  DP4[2]<5 || DP4[3]<5 ||(DP4[2] + DP4[3])/sum(DP4) < 0.9 ||  MQ<50' {input} | "
        "bcftools query  -i'FILTER=\"PASS\"' -f '%LINE' -o {output.fvcf} 2> {log.log}"



rule annotateVars:
    input:
        fvcf = OUTDIR/'VCF/{sample}/{sample}.AF.filtered.vcf',
        #snpEff_gbk =f'{config["snpEff_nenv"]}/data/{config["snpEff_reference"]}/genes.gbk'
    output: avcf = OUTDIR/'VCF/{sample}/{sample}.filtered.annotated.vcf',
        marker = touch(OUTDIR/'VCF/{sample}/{sample}.snpEff.done')
    params:
        genome = config["snpEff_reference"],
        rvcf = lambda wildcards: OUTDIR/f'VCF/{wildcards.sample}/{wildcards.sample}.filtered.renamed.vcf',
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}/{wildcards.sample}.snpEff.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}/{wildcards.sample}.snpEff.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/anVar.yaml'
    log: log = OUTDIR/'logs/VCF/{sample}/{sample}.snpEff.log'
    threads:
        8
    shell:
        'snpEff {params.genome} {input.fvcf} > {output.avcf} 2> {log.log}'


rule buildReferenceSnpEff:
    input: gbk = config["snpEff_gbk"]
    output:
        touch(OUTDIR/f'VCF/{config["snpEff_reference"]}.snpEff_db.done'),
        snpEff_gbk=f'{config["snpEff_nenv"]}/data/{config["snpEff_reference"]}/genes.gbk'
    params:
        genome_name = config["snpEff_reference"],
        fa = config['snpEff_fa'],
        nenv = config['snpEff_nenv'],
        qerrfile = lambda wildcards: OUTDIR/f'logs/{config["snpEff_reference"]}_db.snpEff.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{config["snpEff_reference"]}.snpEff.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/anVar.yaml'
    threads:
        8
    log:
        log = OUTDIR/f'logs/{config["snpEff_reference"]}_db.snpEff.log'

    shell:
        "python scripts/snpEff_db.py {input.gbk} {params.genome_name} "
        " -fa {params.fa} -nenv {params.nenv} ; "
        "snpEff build -genbank -v {params.genome_name}; "
        "cp -r /nfs/home/ansintsova/./miniconda3/pkgs/snpeff-5.0-0/share/snpeff-5.0-0/data {params.nenv}  "




#
# rule annotateVars2:
#     input:
#         fvcf = OUTDIR/'VCF/{sample}/{sample}.AF.filtered.vcf'
#     output: avcf = OUTDIR/'VCF/{sample}/{sample}.AF.filtered.annotated.vcf',
#         marker = touch(OUTDIR/'VCF/{sample}/{sample}.AF.snpEff.done')
#     params:
#         ref = config['reference'],
#         rvcf = lambda wildcards: OUTDIR/f'VCF/{wildcards.sample}/{wildcards.sample}.filtered.renamed.vcf',
#         qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}/{wildcards.sample}.AF.snpEff.qerr',
#         qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}/{wildcards.sample}.AF.snpEff.qout',
#         scratch = 6000,
#         mem = 7700,
#         time = 1400
#     conda:
#         'envs/call_variants.yaml'
#     log: log = OUTDIR/'logs/VCF/{sample}/{sample}.AF.snpEff.log'
#     threads:
#         8
#     shell:
#         'python scripts/snpEff_db.py rename {input.fvcf} {params.rvcf} {params.ref}; '
#         'snpEff {params.ref} {params.rvcf} > {output.avcf} 2> {log.log}'