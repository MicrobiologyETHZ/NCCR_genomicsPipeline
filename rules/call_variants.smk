if config['reference'].endswith('.fasta') or config['reference'].endswith('.fasta.gz'):
    REFPATH = Path(config['reference'])
else:
    REFPATH = OUTDIR/f'assembly/{config["reference"]}/scaffolds.fasta'




rule remove_duplicates:
    input: OUTDIR/'{bam}/{sample}/{sample}.bam'
    output: OUTDIR/'{bam}/{sample}/{sample}.rmdup.bam'
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
    threads:
        8
    shell: "gatk --java-options '-Xmx{params.ram}G' "
           "MarkDuplicates --TMP_DIR {params.tmpdir} "
           "--CREATE_INDEX true --REMOVE_DUPLICATES true "
           "-O {output} -I {input} -M {params.metrics} "



rule bcf_call:
    input:
        scaf = REFPATH,
        bam = OUTDIR/"ref_bams/{sample}/{sample}.bam"
    output:
        vcf = OUTDIR/'VCF/{sample}/{sample}.vcf'
    params:
        qerrfile = OUTDIR/'VCF/{sample}/{sample}.mpileup.qerr',
        qoutfile = OUTDIR/'VCF/{sample}/{sample}.mpileup.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/call_variants.yaml'
    threads:
        8
    shell:
        "bcftools mpileup -Ou -f {input.scaf}  {input.bam} | "
        "bcftools call -mv -o {output.vcf}"


rule bcf_filter:
    input: OUTDIR/'VCF/{sample}/{sample}.vcf'
    output: fvcf = OUTDIR/'VCF/{sample}/{sample}.filtered.vcf',
        marker = touch(OUTDIR/'VCF/{sample}/{sample}.vcf.done')
    params:
        qerrfile = OUTDIR/'VCF/{sample}/{sample}.bcf.qerr',
        qoutfile = OUTDIR/'VCF/{sample}/{sample}.bcf.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/call_variants.yaml'
    threads:
        8
    shell:
        "bcftools filter -Ov -sLowQual -g5 -G10 -e 'QUAL<100 || DP4[2]<10 || DP4[3]<10 ||  MQ<60' {input} | "
        "bcftools query  -i'FILTER=\"PASS\"' -f '%LINE' -o {output.fvcf}"


#DP4[0]>20 || DP4[1]>20 ||

# rule bcf_call:
#          input: OUTDIR/'VCF/{sample}/{sample}.mpileup',
#          output: OUTDIR/'VCF/{sample}/{sample}.vcf'
#          params:
#              qerrfile = OUTDIR/'VCF/{sample}/{sample}.mpileup.qerr',
#              qoutfile = OUTDIR/'VCF/{sample}/{sample}.mpileup.qout',
#              scratch = 6000,
#              mem = 7700,
#              time = 1400
#          conda:
#             'envs/call_variants.yaml'
#          threads:
#               8
#          shell:
#             'bcftools call -m -o {output} {input}'
