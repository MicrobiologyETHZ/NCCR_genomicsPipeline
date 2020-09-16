
rule mlst:
    input:
         scaf = OUTDIR/'assembly/{sample}/scaffolds.fasta'
    output:
         mlst = OUTDIR/'typing/{sample}/{sample}.mlst',
    params:
         scratch = 1000,
         mem = 8000,
         time = 235,
         qerrfile = OUTDIR /'typing/{sample}/{sample}.mlst.qerr',
         qoutfile = OUTDIR /'typing/{sample}/{sample}.mlst.qout'
    conda:
        'envs/typing.yaml'
    shell:
         'mlst {input.scaf} > {output.mlst}'

rule ectyper:
    input:
        scaf = OUTDIR/'assembly/{sample}/scaffolds.fasta'
    output:
        mlst = OUTDIR/'typing/{sample}/output.tsv',
    params:
        outDir = lambda wildcards: OUTDIR/f'typing/{wildcards.sample}',
        scratch = 1000,
        mem = 8000,
        time = 235,
        qerrfile = OUTDIR /'typing/{sample}/{sample}.ectyper.qerr',
        qoutfile = OUTDIR /'typing/{sample}/{sample}.ectyper.qout'
    conda:
        'envs/typing.yaml'
    shell:
         'ectyper -i {input.scaf} -o {params.outDir}'