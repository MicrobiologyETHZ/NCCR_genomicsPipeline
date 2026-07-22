"""
ISMapper — locate insertion sequence (IS) elements relative to a reference genome.

Per sample, maps reads against each IS query and reports the site and orientation
of every insertion, split into "novel" (absent from the reference) and "known"
calls. Runs on the clean reads from `qc` (preprocess.smk) by default, against the
same GenBank reference breseq uses.

Config (see configs/ismap_config.yaml):

    reference:
      refgbk: ../path/to/reference.gbk   # GenBank, shared with breseq
    ismap:
      queries: ../path/to/is_query.fasta # multi-FASTA of IS sequences (required)
      use_raw_reads: false               # true = skip QC, use samples.csv directly
      min_clip: 10       # OPTIONAL, ISMapper defaults are used when omitted
      max_clip: 30
      cutoff: 6
      merging: 100
      all_hits: false    # --a: report all BWA alignments, not just the best

REFGBK is defined in call_variants.smk; OUTDIR, sampleInfo and _resolve in the
Snakefile.

Note on read names: ISMapper derives the sample name from the read filename and
only accepts `<prefix>_1.fastq[.gz]` / `_R1.fastq[.gz]` (see read_grouping.py in
ISMapper). Our clean reads are `{sample}.1.fq.gz`, which matches none of its
patterns, so the rule stages symlinks under the accepted names first. That also
fixes the output directory to `{output_dir}/{sample}/`.

Two logs are written per sample under OUTDIR/logs/ismap/: `{sample}.ismap.log`
(stdout/stderr of the job) and `{sample}.ismapper.log` (ISMapper's own logger).
"""

_ismap_cfg = config.get('ismap', {})
ISMAP_QUERIES = str(_resolve(_ismap_cfg['queries'])) if _ismap_cfg.get('queries') else ''
ISMAP_USE_RAW = _ismap_cfg.get('use_raw_reads', False)

# Flags are only emitted when set, so anything left out of the config falls
# through to ISMapper's own default (same approach as BRESEQ_* in call_variants.smk).
ISMAP_MIN_CLIP = f'--min_clip {_ismap_cfg["min_clip"]}' if _ismap_cfg.get('min_clip') else ''
ISMAP_MAX_CLIP = f'--max_clip {_ismap_cfg["max_clip"]}' if _ismap_cfg.get('max_clip') else ''
ISMAP_CUTOFF = f'--cutoff {_ismap_cfg["cutoff"]}' if _ismap_cfg.get('cutoff') else ''
ISMAP_MERGING = f'--merging {_ismap_cfg["merging"]}' if _ismap_cfg.get('merging') else ''
ISMAP_ALL = '--a' if _ismap_cfg.get('all_hits') else ''
ISMAP_OPTS = ' '.join(o for o in (ISMAP_MIN_CLIP, ISMAP_MAX_CLIP, ISMAP_CUTOFF,
                                  ISMAP_MERGING, ISMAP_ALL) if o)


def _ismap_fq(sample, read):
    """Path to read 1 or 2 for a sample: raw from samples.csv, or QC output.

    getFastq1/getFastq2 in preprocess.smk are not reused here because they are
    only defined when sampleInfo is non-empty.
    """
    if ISMAP_USE_RAW:
        if sampleInfo is None:
            raise WorkflowError(
                "ismap.use_raw_reads needs the read paths from a sample sheet. "
                "Set 'samples: /path/to/samples.csv' in the config, or set "
                "use_raw_reads: false to run on the pipeline's clean reads.")
        column = f'fastq_{read}'
        return str(sampleInfo[sampleInfo['sample'] == sample][column].iloc[0])
    return str(OUTDIR / f'clean_reads/{sample}/{sample}.{read}.fq.gz')


def _ismap_queries(wildcards):
    """The IS query FASTA, validated.

    Checked here rather than at module level so that the error only surfaces for
    runs that actually want ISMapper — this file is parsed for every command.
    """
    if not ISMAP_QUERIES:
        raise WorkflowError(
            "ISMapper needs a multi-FASTA of IS sequences. Set it in the config:\n\n"
            "    ismap:\n      queries: /path/to/is_query.fasta\n")
    if not str(REFGBK) or str(REFGBK) == '.':
        raise WorkflowError(
            "ISMapper needs an annotated GenBank reference. Set it in the config:\n\n"
            "    reference:\n      refgbk: /path/to/reference.gbk\n")
    return ISMAP_QUERIES


rule run_ismap:
    input:
        fq1 = lambda wildcards: _ismap_fq(wildcards.sample, 1),
        fq2 = lambda wildcards: _ismap_fq(wildcards.sample, 2),
        queries = _ismap_queries,
        ref = REFGBK
    output:
        marker = touch(OUTDIR / 'ismap/{sample}/.ismap.done')
    params:
        # ISMapper creates the per-sample subdirectory itself, named after the
        # read prefix, so every sample shares one output_dir.
        outdir = OUTDIR / 'ismap',
        staged = lambda wildcards: OUTDIR / f'ismap/.staged_reads/{wildcards.sample}',
        # Absolute and per-sample: ISMapper's default log name is a timestamp,
        # which concurrent jobs would collide on. ISMapper appends '.log' and
        # opens it with filemode='w', so this prefix must NOT collide with
        # {log.log} below — that file already has the shell redirect writing to it.
        log_prefix = lambda wildcards: OUTDIR / f'logs/ismap/{wildcards.sample}.ismapper',
        opts = ISMAP_OPTS,
        qerrfile = lambda wildcards: OUTDIR / f'logs/ismap/{wildcards.sample}.ismap.qerr',
        qoutfile = lambda wildcards: OUTDIR / f'logs/ismap/{wildcards.sample}.ismap.qout',
        scratch = 6000,
        mem = 7700,
        time = 1400
    conda:
        'ismap'
    log:
        log = OUTDIR / 'logs/ismap/{sample}.ismap.log'
    threads:
        8
    shell:
        "mkdir -p {params.staged} $(dirname {log.log}); "
        "ln -sf $(realpath {input.fq1}) {params.staged}/{wildcards.sample}_1.fastq.gz; "
        "ln -sf $(realpath {input.fq2}) {params.staged}/{wildcards.sample}_2.fastq.gz; "
        "ismap --reads {params.staged}/{wildcards.sample}_1.fastq.gz "
        "{params.staged}/{wildcards.sample}_2.fastq.gz "
        "--queries {input.queries} --reference {input.ref} "
        "--output_dir {params.outdir} --log {params.log_prefix} "
        "--t {threads} {params.opts} &> {log.log}"
