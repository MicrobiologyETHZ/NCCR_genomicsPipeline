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
        OUTDIR/'assembly/{sample1}/{sample1}.scaffolds.min500.fasta.gz'
    output:
        OUTDIR/'assembly/{sample1}/{sample1}.scaffolds.min500.fasta'
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
         ref = OUTDIR/'assembly/{sample1}/{sample1}.scaffolds.min500.fasta',
         query = OUTDIR/'assembly/{sample2}/{sample2}.scaffolds.min500.fasta',
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
        'envs/compare_genomes.yaml'
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


rule runPanX:
    input: get_gbks()
    output:
        tree = OUTDIR/f'panX/data/{config["projectName"]}/vis/strain_tree.nwk'
    params:
        rundir = OUTDIR/f'panX/data/{config["projectName"]}',
        cg = config['cg'],
        projectName = config['projectName'],
        panXpath = config['panXpath'],
        qerrfile = OUTDIR/f'logs/panX/{config["projectName"]}.panX.qerr',
        qoutfile = OUTDIR/f'logs/panX/{config["projectName"]}.panX.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/panX.yaml'
    log:
        log = OUTDIR/f'logs/panX/{config["projectName"]}.panX.log',
    threads:
        32
    shell:
        "mkdir -p {params.rundir}/input_GenBank; cp -f {input} {params.rundir}/input_GenBank; "
        "{params.panXpath}/panX.py -fn {params.rundir} -sl {params.projectName} -cg {params.cg} -t 32 &> {log.log} "


rule setupPhylophlanDB:
    input: get_fnas()
    output: marker = OUTDIR/f'db/phylophlan/{config["species"]}.db.done'
    params:
        outDir = OUTDIR/f'db/phylophlan/',
        sp = config['species'],
        qerrfile = OUTDIR/f'logs/phylophlan/{config["species"]}.db.qerr',
        qoutfile = OUTDIR/f'logs/phylophlan/{config["species"]}.db.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'envs/compare_genomes.yaml'
    log:
        log = OUTDIR/f'logs/phylophlan/{config["species"]}.db.log',
    threads:
        16
    shell:
        "mkdir -p {params.outDir}; "
        "phylophlan_setup_database -g s__{params.sp} "
        "-o {params.outDir} "
        "--verbose &> {log.log} "

rule phylophlan:
    input: fnas = get_fnas(),
         marker = OUTDIR/f'db/phylophlan/{config["species"]}.db.done'
    output: OUTDIR/f'phylophlan/{config["projectName"]}/output_isolates/RAxML_bestTree.input_isolates.tre'
    params:
        runDir = OUTDIR/f'phylophlan/{config["projectName"]}',
        entries = len(get_fnas()),
        sp = config['species'],
        dbFolder = OUTDIR/f'db/phylophlan/',
        qerrfile = OUTDIR/f'logs/phylophlan/{config["species"]}.phylophlan.qerr',
        qoutfile = OUTDIR/f'logs/phylophlan/{config["species"]}.phylophlan.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400,
        phyloConfig = config["phyloConfig"]
    conda:
        'envs/compare_genomes.yaml'
    log:
        log = OUTDIR/f'logs/phylophlan/{config["species"]}.phylophlan.log',
    threads:
        32
    shell:
        "mkdir -p {params.runDir}/input_isolates; "
        "mkdir -p {params.runDir}/output_isolates; "
        "cp -f {input.fnas} {params.runDir}/input_isolates; "
        "phylophlan -i {params.runDir}/input_isolates "
        "-o {params.runDir}/output_isolates "
        "--databases_folder {params.dbFolder} "
        "-d s__{params.sp} "
        "--min_num_entries {params.entries} "
        "--trim greedy "
        "--not_variant_threshold 0.99 "
        "--remove_fragmentary_entries "
        "--fragmentary_threshold 0.67 "
        "-t a " ## doing this on amino acid sequences 
        "-f {params.phyloConfig} "
        "--diversity low "
        "--force_nucleotides "
        "--nproc 32 "
        "--verbose &> {log.log} "


rule ariba:
    input: r1 = OUTDIR/'clean_reads/{sample}/{sample}.1.fq.gz',
         r2 = OUTDIR/'clean_reads/{sample}/{sample}.2.fq.gz'
    output: touch(OUTDIR/'ariba/{sample}.ariba.done')
    params:
        outDir = lambda wildcards: OUTDIR/f'ariba/{wildcards.sample}',
        dbFolder = config["aribaDB"],
        qerrfile = lambda wildcards: OUTDIR/f'logs/ariba/{wildcards.sample}.ariba.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/ariba/{wildcards.sample}.ariba.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400,
    conda:
        'envs/compare_genomes.yaml'
    log:
        log = OUTDIR/'logs/ariba/{sample}.ariba.log',
    threads:
        32
    shell:
        "ariba run {params.dbFolder} "
        "{input.r1} {input.r2} "
        "{params.outDir} &> {log.log} "

rule aribaSummary:
    input: [OUTDIR/f'ariba/{sample}.ariba.done' for sample in SAMPLES]
    output: touch(OUTDIR/'ariba/summary.done')
    params:
        outDir = OUTDIR/'ariba/summary',
        files = [OUTDIR/f'ariba/{sample}/report.tsv' for sample in SAMPLES],
        qerrfile = OUTDIR/f'logs/ariba/summary.qerr',
        qoutfile = OUTDIR/f'logs/ariba/summary.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400,
    conda:
        'envs/compare_genomes.yaml'
    log:
        log = OUTDIR/'logs/ariba/summary.log',
    threads:
        32
    shell:
        "ariba summary {params.outDir} "
        "{params.files}  &> {log.log} "



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
        diamondDB = '/nfs/cds/Databases/DIAMOND/nr',
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





