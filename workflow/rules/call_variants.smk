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

REFDIR = Path(config['reference']['refDir'])
REFGBK =  Path(config['reference']['refgbk']) if 'refgbk' in config['reference'].keys() else ''

rule run_breseq:
    input: fq1=OUTDIR / 'clean_reads/{sample}/{sample}.1.fq.gz',
           fq2=OUTDIR / 'clean_reads/{sample}/{sample}.2.fq.gz',
    output: 
            marker = touch(OUTDIR / 'breseq/{sample}.breseq.done')
    params:
        out_dir = lambda wildcards: OUTDIR/f'breseq/{wildcards.sample}',
        gbk = REFGBK,
        qerrfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}/{wildcards.sample}.breseq.qerr',
        qoutfile=lambda wildcards: OUTDIR / f'logs/{wildcards.sample}/{wildcards.sample}.breseq.qout',
        scratch=6000,
        mem=7700,
        time=1400
    conda:
        'call_variants'
    log:
        log = OUTDIR /'logs/{sample}.breseq.log'
    threads:
        4
    shell:
        """
        breseq -l 120 -j 8 -o {params.out_dir} -r {params.gbk} {input.fq1} {input.fq2}

        """



rule index:
    input: '{file}.{fasta}'
    output: '{file}.{fasta}.bwt',
        marker = touch('{file}.{fasta}.index.done')
    params:
        qerrfile = '{file}.{fasta}.bwa.qerr',
        qoutfile = '{file}.{fasta}.bwa.qout',
        scratch = 6000,
        mem = 7700,
        time=1400
    conda:
        'call_variants'
    threads:
        8
    shell: "bwa index {input}"


rule align_to_ref:
    input:
        fq1=OUTDIR / 'clean_reads/{sample}/{sample}.1.fq.gz',
        fq2=OUTDIR / 'clean_reads/{sample}/{sample}.2.fq.gz',
        fa =REFDIR/'{ref}.fasta',
        index = REFDIR/'{ref}.fasta.bwt',
        m = REFDIR/'{ref}.fasta.index.done'
    output:
        marker = touch(OUTDIR / 'bams/{sample}/{sample}_to_{ref}.bwa.done'),
        bam=OUTDIR / 'bams/{sample}/{sample}_to_{ref}.bam',
        bai=OUTDIR / 'bams/{sample}/{sample}_to_{ref}.bam.bai'
    params:
        qerrfile = lambda wildcards: OUTDIR / 'logs' / f'{wildcards.sample}_{wildcards.ref}.bam.qerr',
        qoutfile = lambda wildcards: OUTDIR / 'logs' / f'{wildcards.sample}_{wildcards.ref}.bam.qout',
        scratch=6000,
        mem=7700,
        time=1400
    log:
        log =  OUTDIR / 'logs' / '{sample}_{ref}.bam.log',
    conda:
        'call_variants'
    threads:
        8
    shell:
        "bwa mem  -t 4 -M {input.fa} "
        "{input.fq1} {input.fq2} | samtools sort --reference {input.fa} "
        "-l 9 -@ 4 -O bam -o {output.bam} -T {output.bam}.tmp; samtools index {output.bam} &> {log.log}"


rule remove_duplicates:
    input: bam = OUTDIR / 'bams/{sample}/{sample}_to_{ref}.bam',
        inmar = REFDIR/'{ref}.fasta.index.done'
    output: bam = OUTDIR/'bams/{sample}/{sample}_to_{ref}.rmdup.bam',
        marker = touch(OUTDIR/'bams/{sample}/{sample}_to_{ref}.rmdup.done')
    params:
        qerrfile = lambda wildcards: OUTDIR/f'bams/{wildcards.sample}/{wildcards.sample}_{wildcards.ref}.markdup.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'bams/{wildcards.sample}/{wildcards.sample}_{wildcards.ref}.markdup.qout',
        scratch = 6000,
        mem = 10000,
        time = 1400,
        ram = config['ram'],
        tmpdir = config['tmpdir'],
        metrics = lambda wildcards: OUTDIR/f'bams/{wildcards.sample}/{wildcards.sample}_{wildcards.ref}.picard.metrics'
    conda:
        'call_variants'
    log:
        log = OUTDIR/'logs/{sample}_{ref}.removeDuplicates.log'
    threads:
        8
    shell: "gatk --java-options '-Xmx{params.ram}G' "
           "MarkDuplicates --TMP_DIR {params.tmpdir} "
           "--CREATE_INDEX true --REMOVE_DUPLICATES true "
           "-O {output.bam} -I {input.bam} -M {params.metrics} &> {log.log} "

# # recalibrate quality?
#

rule bcf_call:
    input:
        fa = REFDIR/'{ref}.fasta',
        bam = OUTDIR/"bams/{sample}/{sample}_to_{ref}.rmdup.bam"
    output:
        vcf = OUTDIR/'VCF/{sample}/{sample}_to_{ref}.vcf',
        marker = touch(OUTDIR/'VCF/{sample}/{sample}_to_{ref}.mpileup.done')
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}_{wildcards.ref}.mpileup.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}_{wildcards.ref}.mpileup.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    log:
        log = OUTDIR/'logs/VCF/{sample}_{ref}.mpileup.log'
    conda:
        'call_variants'
    threads:
        32
    shell:
        "bcftools mpileup -Ou -q30 -d3000 -f {input.fa}  {input.bam} | "
        "bcftools call -mv -o {output.vcf} 2> {log.log}"

rule bcf_filter_isolate:
    input: OUTDIR / 'VCF/{sample}/{sample}_to_{ref}.vcf',
    output: fvcf = OUTDIR / 'VCF/{sample}/{sample}_to_{ref}.isolate.filtered.vcf',
            marker = touch(OUTDIR / 'VCF/{sample}/{sample}_to_{ref}.isolate.vcf.done')
    params:
        qerrfile=lambda wildcards: OUTDIR / f'logs/VCF/{wildcards.sample}/{wildcards.sample}_{wildcards.ref}.isolate.bcf.qerr',
        qoutfile=lambda wildcards: OUTDIR / f'logs/VCF/{wildcards.sample}/{wildcards.sample}_{wildcards.ref}.isolate.bcf.qout',
        scratch=6000,
        mem=7700,
        time=1400
    conda:
        'call_variants'
    log:
        log = OUTDIR / 'logs/VCF/{sample}_{ref}.isolate.bcf.log'
    threads:
        8
    shell:
        #"bcftools filter -Ov -sLowQual -g5 -G10 -e 'QUAL<10 ||  DP4[2]<10 || DP4[3]<10 ||(DP4[2] + DP4[3])/sum(DP4) < 0.9 ||  MQ<50' {input} | "
        "bcftools filter -Ov -sLowQual -g5 -G10 -e 'QUAL<10 ||  DP4[2]<5 || DP4[3]<5 ||(DP4[2] + DP4[3])/sum(DP4) < 0.5 ||  MQ<50' {input} | "
        "bcftools query  -i'FILTER=\"PASS\"' -f '%LINE' -o {output.fvcf} &> {log.log}"


rule annotateVars:
    input:
        fvcf = OUTDIR/'VCF/{sample}/{sample}.AF.filtered.vcf',
    output:
        avcf = OUTDIR/'VCF/{sample}/{sample}.filtered.annotated.vcf',
        marker = touch(OUTDIR/'VCF/{sample}/{sample}.snpEff.done')
    params:
        genome = config["snpEff_reference"],
        #rvcf = lambda wildcards: OUTDIR/f'VCF/{wildcards.sample}/{wildcards.sample}.filtered.renamed.vcf',
        qerrfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}/{wildcards.sample}.snpEff.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/{wildcards.sample}/{wildcards.sample}.snpEff.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'anVar'
    log: log = OUTDIR/'logs/VCF/{sample}/{sample}.snpEff.log'
    threads:
        8
    shell:
        'snpEff {params.genome} {input.fvcf} > {output.avcf} 2> {log.log}'


# rule bcf_filter:
#     input: vcf = OUTDIR/'VCF/{sample}/{sample}.vcf',
#         m = OUTDIR/'VCF/{sample}/{sample}.mpileup.done'
#     output: fvcf = OUTDIR/'VCF/{sample}/{sample}.filtered.vcf',
#         marker = touch(OUTDIR/'VCF/{sample}/{sample}.vcf.done')
#     params:
#         qerrfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}.bcf.qerr',
#         qoutfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}.bcf.qout',
#         scratch = 6000,
#         mem = 7700,
#         time = 1400
#     log:
#         log = OUTDIR/'logs/VCF/{sample}.bcf.log'
#     conda:
#         'envs/call_variants.yaml'
#     threads:
#         8
#     shell:
#         "bcftools filter -Ov -sLowQual -g5 -G10 -e 'QUAL<100 || DP4[2]<10 || DP4[3]<10 ||  MQ<60' {input.vcf} | "
#         "bcftools query  -i'FILTER=\"PASS\"' -f '%LINE' -o {output.fvcf} 2> {log.log}"
#
#

#
#
# rule bcf_filter3:
#     input: OUTDIR/'VCF/{sample}/{sample}.vcf',
#
#     output: fvcf = OUTDIR/'VCF/{sample}/{sample}.f3.filtered.vcf',
#         marker = touch(OUTDIR/'VCF/{sample}/{sample}.f3.vcf.done')
#     params:
#         qerrfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}.f3.bcf.qerr',
#         qoutfile = lambda wildcards: OUTDIR/f'logs/VCF/{wildcards.sample}/{wildcards.sample}.f3.bcf.qout',
#         scratch = 6000,
#         mem = 7700,
#         time = 1400
#     conda:
#         'envs/call_variants.yaml'
#     log:
#         log = OUTDIR/'logs/VCF/{sample}.f3.bcf.log'
#     threads:
#         8
#     shell:
#         "bcftools filter -Ov -sLowQual -g5 -G10 -e 'QUAL<10 ||  DP4[2]<5 || DP4[3]<5 ||(DP4[2] + DP4[3])/sum(DP4) < 0.9 ||  MQ<50' {input} | "
#         "bcftools query  -i'FILTER=\"PASS\"' -f '%LINE' -o {output.fvcf} 2> {log.log}"
#
#
#

#
#
# rule buildReferenceSnpEff:
#     input: gbk = config["snpEff_gbk"]
#     output:
#         touch(OUTDIR/f'VCF/{config["snpEff_reference"]}.snpEff_db.done'),
#         snpEff_gbk=f'{config["snpEff_nenv"]}/data/{config["snpEff_reference"]}/genes.gbk'
#     params:
#         genome_name = config["snpEff_reference"],
#         fa = config['snpEff_fa'],
#         nenv = config['snpEff_nenv'],
#         qerrfile = lambda wildcards: OUTDIR/f'logs/{config["snpEff_reference"]}_db.snpEff.qerr',
#         qoutfile = lambda wildcards: OUTDIR/f'logs/{config["snpEff_reference"]}.snpEff.qout',
#         scratch = 6000,
#         mem = 7700,
#         time = 1400
#     conda:
#         'envs/anVar.yaml'
#     threads:
#         8
#     log:
#         log = OUTDIR/f'logs/{config["snpEff_reference"]}_db.snpEff.log'
#
#     shell:
#         "python scripts/snpEff_db.py {input.gbk} {params.genome_name} "
#         " -fa {params.fa} -nenv {params.nenv} ; "
#         "snpEff build -genbank -v {params.genome_name}; "
#         "cp -r /nfs/home/ansintsova/./miniconda3/pkgs/snpeff-5.0-0/share/snpeff-5.0-0/data {params.nenv}  "
#
#


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