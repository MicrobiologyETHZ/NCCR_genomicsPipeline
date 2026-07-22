"""
InStrain — strain-level microdiversity profiling and cross-sample comparison.

Consumes the sorted/indexed BAMs produced by `align_to_ref` (call_variants.smk),
i.e. reads from each sample mapped to a reference FASTA in REFDIR. Per (sample, ref)
it runs `inStrain profile`; per reference with >=2 samples it runs `inStrain compare`
to get population-level ANI / strain sharing.

Config (see configs/instrain_config.yaml):

    reference:
      refDir: ../path/to/refs     # directory holding {ref}.fasta
      myref: [sample1, sample2]   # reference name -> samples mapped to it
    instrain:
      min_cov: 5                  # --min_cov for profiling
      run_compare: true           # also run inStrain compare across samples
      genes:                      # OPTIONAL, keyed by reference name
        myref: ../path/to/myref.genes.fna   # prodigal-style gene calls
      stb:                        # OPTIONAL, keyed by reference name
        myref: ../path/to/myref.stb         # scaffold-to-bin (genome-level metrics)

REFDIR is defined in call_variants.smk; _resolve is defined in the Snakefile.
"""

_instrain_cfg = config.get('instrain', {})
IS_MIN_COV = _instrain_cfg.get('min_cov', 5)
# Resolve user-supplied gene/stb files (relative paths are relative to workflow/).
IS_GENES = {ref: str(_resolve(p)) for ref, p in _instrain_cfg.get('genes', {}).items()}
IS_STB = {ref: str(_resolve(p)) for ref, p in _instrain_cfg.get('stb', {}).items()}


def _genes_input(wildcards):
    """Optional gene-calls file for this reference, or [] if none provided."""
    g = IS_GENES.get(wildcards.ref)
    return [g] if g else []


def _stb_input(wildcards):
    """Optional scaffold-to-bin file for this reference, or [] if none provided."""
    s = IS_STB.get(wildcards.ref)
    return [s] if s else []


def _compare_profile_markers(wildcards):
    """Profile-done markers for every sample mapped to this reference."""
    samples = config['reference'][wildcards.ref]
    return [OUTDIR / f'instrain/profiles/{s}_to_{wildcards.ref}/.IS.profile.done'
            for s in samples]


rule run_instrain_profile:
    input:
        bam = OUTDIR / 'bams/{sample}/{sample}_to_{ref}.bam',
        bai = OUTDIR / 'bams/{sample}/{sample}_to_{ref}.bam.bai',
        fa = REFDIR / '{ref}.fasta',
        genes = _genes_input,
        stb = _stb_input,
    output:
        marker = touch(OUTDIR / 'instrain/profiles/{sample}_to_{ref}/.IS.profile.done')
    params:
        outdir = lambda wildcards: OUTDIR / f'instrain/profiles/{wildcards.sample}_to_{wildcards.ref}',
        genes_arg = lambda wildcards, input: f'-g {input.genes[0]}' if input.genes else '',
        stb_arg = lambda wildcards, input: f'-s {input.stb[0]}' if input.stb else '',
        min_cov = IS_MIN_COV,
        qerrfile = lambda wildcards: OUTDIR / f'logs/instrain/{wildcards.sample}_to_{wildcards.ref}.profile.qerr',
        qoutfile = lambda wildcards: OUTDIR / f'logs/instrain/{wildcards.sample}_to_{wildcards.ref}.profile.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'instrain'
    log:
        log = OUTDIR / 'logs/instrain/{sample}_to_{ref}.profile.log'
    threads:
        8
    shell:
        "inStrain profile {input.bam} {input.fa} -o {params.outdir} "
        "-p {threads} {params.genes_arg} {params.stb_arg} --min_cov {params.min_cov} &> {log.log}"


rule run_instrain_compare:
    input:
        markers = _compare_profile_markers,
        stb = _stb_input,
    output:
        marker = touch(OUTDIR / 'instrain/compare/{ref}/.IS.compare.done')
    params:
        profiles = lambda wildcards: ' '.join(
            str(OUTDIR / f'instrain/profiles/{s}_to_{wildcards.ref}')
            for s in config['reference'][wildcards.ref]),
        stb_arg = lambda wildcards, input: f'-s {input.stb[0]}' if input.stb else '',
        outdir = lambda wildcards: OUTDIR / f'instrain/compare/{wildcards.ref}',
        qerrfile = lambda wildcards: OUTDIR / f'logs/instrain/{wildcards.ref}.compare.qerr',
        qoutfile = lambda wildcards: OUTDIR / f'logs/instrain/{wildcards.ref}.compare.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'instrain'
    log:
        log = OUTDIR / 'logs/instrain/{ref}.compare.log'
    threads:
        8
    shell:
        "inStrain compare -i {params.profiles} -o {params.outdir} "
        "-p {threads} {params.stb_arg} &> {log.log}"
