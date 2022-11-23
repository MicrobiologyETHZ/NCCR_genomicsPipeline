


rule motus_profile:
    input:
        fq1 = OUTDIR/READ_DIR/'{sample}/{sample}.1.fq.gz',
        fq2 = OUTDIR/READ_DIR/'{sample}/{sample}.2.fq.gz',
        s = OUTDIR/READ_DIR/'{sample}/{sample}.s.fq.gz',
    output:
        OUTDIR/'motus/{sample}/{sample}.motus'
    params:
          qerrfile=lambda wildcards: OUTDIR / f'logs/motus/{wildcards.sample}.qerr',
          qoutfile=lambda wildcards: OUTDIR / f'logs/motus/{wildcards.sample}.qout',
          scratch=6000,
          mem=7700,
          time=1400
    log:
       log=OUTDIR/'logs/motus/{sample}.log',
    conda:
         'envs/profile.yaml'
    threads:
           8
    shell:
         "motus profile -f {input.fq1} -r {input.fq2} -s {input.s} -p -c -o {output}"


rule fetchMG_profile:
    input:
        fna = OUTDIR/'assembly/{sample}/prokka/{sample}.fna',
        faa = OUTDIR/'assembly/{sample}/prokka/{sample}.faa'
    output:
        #OUTDIR/'fetchMG/{sample}/marker_genes_scores.table'
        touch(OUTDIR/'fetchMG/{sample}/fetchMG.done')
    params:
          outdir = lambda wildcards: OUTDIR / f'fetchMG/{wildcards.sample}/output',
          fetchMGpath = config['pathfetchMG'],
          qerrfile=lambda wildcards: OUTDIR / f'logs/fetchMG/{wildcards.sample}.qerr',
          qoutfile=lambda wildcards: OUTDIR / f'logs/fetchMG/{wildcards.sample}.qout',
          scratch=6000,
          mem=7700,
          time=1400
    log:
       log=OUTDIR/'logs/fetchMG/{sample}.log',
    conda:
         'envs/profile.yaml'
    threads:
           8
    shell:
         "{params.fetchMGpath}/fetchMGs.pl -m extraction "
         "-x {params.fetchMGpath}/bin {input.faa} -d {input.fna} -o {params.outdir} "



rule mash:
    input:
        fq1 = OUTDIR/READ_DIR/'{sample}/{sample}.1.fq.gz',
        fq2 = OUTDIR/READ_DIR/'{sample}/{sample}.2.fq.gz',
    output:
        zfq = temp(OUTDIR/'mash/{sample}/{sample}.cat.fastq'),
        msh = temp(OUTDIR/'mash/{sample}/{sample}.cat.fastq.msh'),
        dist = OUTDIR/'mash/{sample}/{sample}.distances.tab',

    params:
          refseqMSH = config['refseqMSH'],
          qerrfile=lambda wildcards: OUTDIR / f'logs/mash/{wildcards.sample}.qerr',
          qoutfile=lambda wildcards: OUTDIR / f'logs/mash/{wildcards.sample}.qout',
          scratch=6000,
          mem=7700,
          time=1400
    log:
       log=OUTDIR/'logs/mash/{sample}.log',
    conda:
         'envs/profile.yaml'
    threads:
           8
    shell:
         "zcat {input.fq1} {input.fq2} > {output.zfq}; "
         "mash sketch -m 2 {output.zfq}; "
         "mash dist {params.refseqMSH} {output.msh} > {output.dist} 2> {log.log}"

