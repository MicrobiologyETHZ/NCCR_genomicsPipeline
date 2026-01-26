from pathlib import Path

# DATADIR = Path(config['dataDir'])
# OUTDIR = Path(config['outDir'])
# REFSTRAIN = OUTDIR/'assembly/LL6/LL6.scaffolds.min500.fasta'
# SFILE = Path(config['sampleFile'])
# SUBSAMPLES = gv.get_subsamples(SFILE)
# SAMPLES = gv.get_samples(DATADIR, SUBSAMPLES)
# rule align_genomes:
#     input: [OUTDIR/f'mummer/{subsample}/{subsample}.report' for subsample in SUBSAMPLES]

# rule dnadiff:
#     input:
#          ref = OUTDIR/'assembly/{sample1}/{sample1}.scaffolds.min500.fasta',
#          query = OUTDIR/'assembly/{sample2}/{sample2}.scaffolds.min500.fasta',
#     output:
#         OUTDIR/'mummer/{sample}/{sample}.report'
#     params:
#         sample = '{sample}',
#         prefix = OUTDIR/'mummer/',
#         qerrfile = OUTDIR/'mummer/{sample}.dnadiff.qerr',
#         qoutfile = OUTDIR/'mummer/{sample}.dnadiff.qout',
#         scratch = 6000,
#         mem = 7700,
#         time = 1400
#     conda:
#         'envs/compare_genomes.yaml'
#     threads:
#         8
#     shell:
#         "dnadiff -p {params.prefix}/{params.sample}/{params.sample} {input.ref} {input.query}"


rule gunzipNucmer:
    input:
        OUTDIR/'assembly/{sample1}/{sample1}.scaffolds.min200.fasta.gz'
    output:
        OUTDIR/'assembly/{sample1}/{sample1}.scaffolds.min200.fasta'
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/mummer/{wildcards.sample1}.gzip.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/mummer/{wildcards.sample1}.gzip.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/compare_genomes.yaml'
    threads:
        8
    shell:
        "gunzip {input} "


rule run_nucmer:
    input:
         ref = OUTDIR/'assembly/{sample1}/{sample1}.scaffolds.min200.fasta',
         query = OUTDIR/'assembly/{sample2}/{sample2}.scaffolds.min200.fasta',
    output:
        OUTDIR/'mummer/{sample1}_{sample2}.delta'
    params:
        prefix = lambda wildcards: OUTDIR/f'mummer/{wildcards.sample1}_{wildcards.sample2}',
        qerrfile = lambda wildcards: OUTDIR/f'logs/mummer/{wildcards.sample1}_{wildcards.sample2}.nucmer.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/mummer/{wildcards.sample1}_{wildcards.sample2}.nucmer.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/compare_genomes.yaml'
    threads:
        8
    log:
        log = OUTDIR/'logs/mummer/{sample1}_{sample2}.nucmer.log'
    shell:
        "nucmer -c 100 -p {params.prefix} {input.ref} {input.query} &> {log.log}"


rule showCoords:
    input:  OUTDIR/'mummer/{sample1}_{sample2}.delta'
    output: OUTDIR/'mummer/{sample1}_{sample2}.coords'
    params:
        prefix = lambda wildcards: OUTDIR/f'mummer/{wildcards.sample1}_{wildcards.sample2}',
        qerrfile = lambda wildcards: OUTDIR/f'logs/mummer/{wildcards.sample1}_{wildcards.sample2}.showcoords.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/mummer/{wildcards.sample1}_{wildcards.sample2}.showcoords.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/compare_genomes.yaml'
    threads:
        8
    log:
        log = OUTDIR/'logs/mummer/{sample1}_{sample2}.nucmerCoords.log'
    shell:
        "show-coords -r -l -c  {input} > {output} 2> {log.log}"



rule calculateANI:
    input: ref = OUTDIR/'{assembly}/{sample1}/scaffolds.fasta.gz',
         genome = OUTDIR/'{assembly}/{sample2}/scaffolds.fasta.gz'
    output:
        aniFile = OUTDIR/'ANI/{assembly}/{sample1}_{sample2}_fastani.out'
    params:
        qerrfile = lambda wildcards: OUTDIR/f'logs/ANI/{wildcards.sample1}_{wildcards.sample2}.fastani.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/ANI/{wildcards.sample1}_{wildcards.sample2}.fastani.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        compgen
    threads:
        8
    log:
        log = OUTDIR/'logs/ANI/{assembly}/{sample1}_{sample2}.fastani.log'
    shell:
        "fastANI -q {input.genome} -r {input.ref} --visualize -o {output.aniFile} &> {log.log}"


def get_gbks():
    if config["gbkFolder"]:
        gbkPath = Path(config["gbkFolder"])
    else:
        gbkPath = Path(OUTDIR/"assembly")
    return [gbk for gbk in gbkPath.glob('**/*.gbk')]

def get_fnas():
    if config['fnaFolder']:
        fnaFolder = Path(config['fnaFolder'])
    else:
        fnaFolder = Path(OUTDIR/"assembly")
    return [fna for fna in fnaFolder.glob('**/*.fna*')]


rule iqtree:
    input: '{msa_file}'
    output: '{msa_file}.iqtree'
    params:
        qerrfile = lambda wildcards: f'{wildcards.msa_file}.qerr',
        qoutfile = lambda wildcards: f'{wildcards.msa_file}.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400,
    conda:
        'envs/compare_genomes.yaml'
    log:
        log = '{msa_file}.log'
    threads:
        32
    shell:
       # "iqtree -s {input} -mset WAG,LG,DCmut -T AUTO -B 1000"
        "iqtree -s {input} -m LG+F+R  -T AUTO -B 1000"


rule calculat_genlen:
    input: '{assembly}'
    output: touch('{assembly}.genelen.done')
    params:
        qerrfile = '{assembly}.genelen.done.qerr',
        qoutfile = '{assembly}.genelen.done.qout',
        diamondDB = config['diamondDB'],
        scratch = 6000,
        mem = 7700,
        time = 1400,
    conda:
        'envs/annotate.yaml'
    log:
        log = '{assembly}.log'
    threads:
        16
    shell:
        "python ./scripts/genlen.py {input} --blast --db {params.diamondDB} "
#
# rule orthoFinder:
#     input: config['orthoFinderData']
#     output: touch(OUTDIR/'orthoFinder.done')
#     params:
#         qerrfile = OUTDIR/'logs/orthofinder.qerr',
#         qoutfile = OUTDIR/'logs/orthofinder.qout',
#         scratch = 6000,
#         mem = 7700,
#         time = 1400,
#     conda:
#         'envs/orthoFinder.yaml'
#     log:
#         log = OUTDIR/'logs/orthofinder.log'
#     threads:
#         16
#     shell:
#         "orthofinder -f {input}"


#
# rule gunzip_assembly:
#     input: OUTDIR/'assembly/{sample}.fasta.gz'
#     output: OUTDIR/'assembly/{sample}.fasta'
#     params:
#         qerrfile = OUTDIR/'assembly/{sample}.gzip.qerr',
#         qoutfile = OUTDIR/'assembly/{sample}.gzip.qout',
#         scratch = 6000,
#         mem = 7700,
#         time = 1400
#     conda:
#         'envs/mummer.yaml'
#     threads:
#         8
#     shell:
#         'gunzip {input}'





