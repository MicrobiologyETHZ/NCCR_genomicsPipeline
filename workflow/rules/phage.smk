"""
Phage / provirus prediction and annotation.

Prediction runs two independent callers — geNomad and Cenote-Taker 3 — and each
caller's viral contigs are then annotated through the same three-step chain,
pharokka -> phold -> phynteny_transformer, so the two callers can be compared
directly. CheckV runs alongside for completeness/contamination and feeds the
summary tables; it does not filter contigs before annotation.

Everything is keyed on two wildcards: `{name}` (an assembly, from
scripts/phage_inputs.resolve_assemblies) and `{caller}` (genomad|cenotetaker).
Outputs live under OUTDIR/phage/, never next to the input assembly, so the
workflow works on read-only collaborator or public data.

Expects these names from the including Snakefile (see Snakefile_phage):
    OUTDIR, ASSEMBLIES, CALLERS, DBS, is_compressed
"""
from pathlib import Path

# The repo has no other wildcard_constraints, and several existing rules use
# greedy free wildcards (compare_genomes.smk:135,153). Constraining ours keeps
# the phage rules from being drawn into those ambiguities.
wildcard_constraints:
    name = r'[A-Za-z0-9._\-]+',
    caller = r'genomad|cenotetaker'


_genomad_cfg = config.get('genomad', {})
_ct_cfg = config.get('cenotetaker', {})
_pharokka_cfg = config.get('pharokka', {})
_phold_cfg = config.get('phold', {})

# Every boolean goes through as_bool: `--config disable_nn=false` arrives as the
# string "false", which is truthy, so a plain `.get()` would silently invert it.
GENOMAD_NN_ARG = '--disable-nn-classification' if as_bool(_genomad_cfg.get('disable_nn'), True) else ''
CT_PROPHAGE = 'T' if as_bool(_ct_cfg.get('prophage'), True) else 'F'
CT_MIN_HALLMARK = _ct_cfg.get('min_hallmark_genes', 1)
PHAROKKA_META = '-m' if as_bool(_pharokka_cfg.get('meta'), True) else ''
PHOLD_CPU = '--cpu' if as_bool(_phold_cfg.get('cpu'), True) else ''

# Decompressed/normalised inputs are intermediates; keep them only if asked.
_intermediate = (lambda p: p) if as_bool(config.get('keep_decompressed'), False) else temp


rule phage_prepare_input:
    """Normalise every assembly to OUTDIR/phage/{name}/input/{name}.fasta.

    Both gzipped and plain inputs are materialised here rather than read in
    place. That costs one genome-sized file per assembly, but it makes every
    downstream tool's output filenames deterministic: geNomad, Cenote-Taker and
    pharokka all derive their output prefixes from the input basename, so
    feeding them `{name}.fasta` is what lets the rules below declare real output
    paths instead of only a marker file.

    Note this decompresses to a new file rather than calling `gunzip` in place
    like alignment.smk:69-81 does — that rule deletes its own source, which
    would be destructive on collaborator data.
    """
    input:
        assembly = lambda wildcards: str(ASSEMBLIES[wildcards.name])
    output:
        fasta = _intermediate(OUTDIR/'phage/{name}/input/{name}.fasta')
    params:
        # Whether the source needs decompressing is known at DAG time, so pick
        # the reader here rather than branching in shell.
        reader = lambda wildcards: ('gzip -cd' if is_compressed(ASSEMBLIES[wildcards.name])
                                    else 'cat'),
        qerrfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.prepare.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.prepare.qout',
        scratch = 1000,
        mem = 2000,
        time = 20
    log:
        log = OUTDIR/'logs/phage/{name}.prepare.log'
    threads:
        1
    shell:
        '{params.reader} {input.assembly} > {output.fasta} 2> {log.log}'


rule genomad:
    """geNomad end-to-end virus/plasmid identification."""
    input:
        fasta = OUTDIR/'phage/{name}/input/{name}.fasta'
    output:
        summary = OUTDIR/'phage/{name}/genomad/{name}_summary/{name}_virus_summary.tsv',
        marker = touch(OUTDIR/'phage/{name}/genomad/{name}.genomad.done')
    params:
        outdir = lambda wildcards: OUTDIR/f'phage/{wildcards.name}/genomad',
        db = DBS['genomad'],
        nn_arg = GENOMAD_NN_ARG,
        qerrfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.genomad.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.genomad.qout',
        scratch = 6000,
        mem = 4000,
        time = 1400
    conda:
        conda_env('genomad')
    log:
        log = OUTDIR/'logs/phage/{name}.genomad.log'
    threads:
        16
    shell:
        # No --cleanup: it removes the intermediate module directories, and
        # collect_viral reads _find_proviruses/ as a fallback source.
        'genomad end-to-end --threads {threads} {params.nn_arg} '
        '{input.fasta} {params.outdir} {params.db} &> {log.log}; '
        # geNomad omits summary files when nothing is classified; guarantee the
        # declared output exists so a virus-free assembly doesn't fail the DAG.
        'mkdir -p $(dirname {output.summary}); touch {output.summary}'


rule cenotetaker:
    """Cenote-Taker 3 virus discovery.

    Cenote-Taker writes its run directory relative to the working directory and
    has no output-directory flag, so the rule cds into its own scratch area
    first. `{input.fasta}` is absolute (it lives under OUTDIR), so the cd is safe.
    """
    input:
        fasta = OUTDIR/'phage/{name}/input/{name}.fasta'
    output:
        summary = OUTDIR/'phage/{name}/cenotetaker/{name}/{name}_virus_summary.tsv',
        marker = touch(OUTDIR/'phage/{name}/cenotetaker/{name}.cenotetaker.done')
    params:
        workdir = lambda wildcards: OUTDIR/f'phage/{wildcards.name}/cenotetaker',
        db = DBS['cenotetaker'],
        prophage = CT_PROPHAGE,
        min_hallmark = CT_MIN_HALLMARK,
        qerrfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.cenotetaker.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.cenotetaker.qout',
        scratch = 6000,
        mem = 4000,
        time = 1400
    conda:
        conda_env('cenotetaker')
    log:
        log = OUTDIR/'logs/phage/{name}.cenotetaker.log'
    threads:
        16
    shell:
        'export CENOTE_DBS={params.db}; '
        'mkdir -p {params.workdir}; cd {params.workdir}; '
        'cenotetaker3 -c {input.fasta} -r {wildcards.name} -p {params.prophage} '
        '-t {threads} --lin_minimum_hallmark_genes {params.min_hallmark} '
        '&> {log.log}; '
        'mkdir -p $(dirname {output.summary}); touch {output.summary}'


rule collect_viral:
    """Normalise each caller's viral contigs to one FASTA with unique headers.

    Contig names are prefixed with the assembly name so they stay unique once
    all assemblies are merged into the cross-assembly summary tables.
    """
    input:
        marker = lambda wildcards: (
            OUTDIR/f'phage/{wildcards.name}/genomad/{wildcards.name}.genomad.done'
            if wildcards.caller == 'genomad'
            else OUTDIR/f'phage/{wildcards.name}/cenotetaker/{wildcards.name}.cenotetaker.done')
    output:
        fasta = OUTDIR/'phage/{name}/viral/{name}.{caller}.fna'
    params:
        # geNomad's summary FASTA is documented to already contain proviruses,
        # but the provirus FASTA is passed too and the script deduplicates by
        # header — so this stays correct either way.
        sources = lambda wildcards: (
            [OUTDIR/f'phage/{wildcards.name}/genomad/{wildcards.name}_summary/{wildcards.name}_virus.fna',
             OUTDIR/f'phage/{wildcards.name}/genomad/{wildcards.name}_find_proviruses/{wildcards.name}_provirus.fna']
            if wildcards.caller == 'genomad'
            else [OUTDIR/f'phage/{wildcards.name}/cenotetaker/{wildcards.name}/{wildcards.name}_virus_sequences.fna']),
        qerrfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.{wildcards.caller}.collect.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.{wildcards.caller}.collect.qout',
        scratch = 1000,
        mem = 2000,
        time = 20
    log:
        log = OUTDIR/'logs/phage/{name}.{caller}.collect.log'
    threads:
        1
    shell:
        'python ./scripts/phage_collect.py --name {wildcards.name} '
        '--output {output.fasta} {params.sources} &> {log.log}'


rule checkv:
    """CheckV completeness/contamination on the predicted viral contigs."""
    input:
        fasta = OUTDIR/'phage/{name}/viral/{name}.{caller}.fna'
    output:
        summary = OUTDIR/'phage/{name}/checkv/{caller}/quality_summary.tsv'
    params:
        outdir = lambda wildcards: OUTDIR/f'phage/{wildcards.name}/checkv/{wildcards.caller}',
        db = DBS['checkv'],
        qerrfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.{wildcards.caller}.checkv.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.{wildcards.caller}.checkv.qout',
        scratch = 6000,
        mem = 4000,
        time = 1400
    conda:
        conda_env('checkv')
    log:
        log = OUTDIR/'logs/phage/{name}.{caller}.checkv.log'
    threads:
        8
    shell:
        'if [ -s {input.fasta} ]; then '
        '  checkv end_to_end {input.fasta} {params.outdir} -t {threads} '
        '  -d {params.db} &> {log.log}; '
        'else '
        '  echo "No viral contigs for {wildcards.name}/{wildcards.caller}; skipping CheckV" > {log.log}; '
        'fi; '
        'mkdir -p {params.outdir}; touch {output.summary}'


rule pharokka:
    """pharokka — step 1 of the annotation chain."""
    input:
        fasta = OUTDIR/'phage/{name}/viral/{name}.{caller}.fna'
    output:
        gbk = OUTDIR/'phage/{name}/annotate/{caller}/pharokka/pharokka.gbk'
    params:
        outdir = lambda wildcards: OUTDIR/f'phage/{wildcards.name}/annotate/{wildcards.caller}/pharokka',
        db = DBS['pharokka'],
        meta = PHAROKKA_META,
        qerrfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.{wildcards.caller}.pharokka.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.{wildcards.caller}.pharokka.qout',
        scratch = 6000,
        mem = 4000,
        time = 1400
    conda:
        conda_env('pharokka')
    log:
        log = OUTDIR/'logs/phage/{name}.{caller}.pharokka.log'
    threads:
        8
    shell:
        'if [ -s {input.fasta} ]; then '
        '  pharokka.py -i {input.fasta} -o {params.outdir} -d {params.db} '
        '  -t {threads} -p pharokka {params.meta} -f &> {log.log}; '
        'else '
        '  echo "No viral contigs for {wildcards.name}/{wildcards.caller}; skipping pharokka" > {log.log}; '
        'fi; '
        'mkdir -p {params.outdir}; touch {output.gbk}'


rule phold:
    """phold — step 2 of the annotation chain, consumes pharokka's GenBank.

    Runs on CPU (`--cpu`); phold's ProstT5 step is much slower without a GPU, so
    this rule gets a long time limit.
    """
    input:
        gbk = OUTDIR/'phage/{name}/annotate/{caller}/pharokka/pharokka.gbk'
    output:
        gbk = OUTDIR/'phage/{name}/annotate/{caller}/phold/phold.gbk'
    params:
        outdir = lambda wildcards: OUTDIR/f'phage/{wildcards.name}/annotate/{wildcards.caller}/phold',
        db = DBS['phold'],
        cpu = PHOLD_CPU,
        qerrfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.{wildcards.caller}.phold.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.{wildcards.caller}.phold.qout',
        scratch = 8000,
        mem = 8000,
        time = 2880
    conda:
        conda_env('phold')
    log:
        log = OUTDIR/'logs/phage/{name}.{caller}.phold.log'
    threads:
        16
    shell:
        'if [ -s {input.gbk} ]; then '
        '  phold run -i {input.gbk} -o {params.outdir} -d {params.db} '
        '  -t {threads} -p phold {params.cpu} -f &> {log.log}; '
        'else '
        '  echo "No pharokka annotation for {wildcards.name}/{wildcards.caller}; skipping phold" > {log.log}; '
        'fi; '
        'mkdir -p {params.outdir}; touch {output.gbk}'


rule phynteny:
    """phynteny_transformer — step 3, refines hypothetical proteins by synteny."""
    input:
        gbk = OUTDIR/'phage/{name}/annotate/{caller}/phold/phold.gbk'
    output:
        gbk = OUTDIR/'phage/{name}/annotate/{caller}/phynteny/phynteny_transformer.gbk',
        cds = OUTDIR/'phage/{name}/annotate/{caller}/phynteny/phynteny_per_cds_funcions.tsv'
    params:
        outdir = lambda wildcards: OUTDIR/f'phage/{wildcards.name}/annotate/{wildcards.caller}/phynteny',
        models = DBS['phynteny'],
        qerrfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.{wildcards.caller}.phynteny.qerr',
        qoutfile = lambda wildcards: OUTDIR/f'logs/phage/{wildcards.name}.{wildcards.caller}.phynteny.qout',
        scratch = 6000,
        mem = 4000,
        time = 1400
    conda:
        conda_env('phynteny')
    log:
        log = OUTDIR/'logs/phage/{name}.{caller}.phynteny.log'
    threads:
        8
    shell:
        # phynteny refuses to write into an existing directory, and Snakemake
        # has already created it for the declared outputs — clear it first.
        'if [ -s {input.gbk} ]; then '
        '  rm -rf {params.outdir}; '
        '  phynteny_transformer {input.gbk} -o {params.outdir} '
        '  -m {params.models} &> {log.log}; '
        'else '
        '  echo "No phold annotation for {wildcards.name}/{wildcards.caller}; skipping phynteny" > {log.log}; '
        'fi; '
        'mkdir -p {params.outdir}; touch {output.gbk} {output.cds}'


rule phage_summary:
    """Merge per-assembly, per-caller results into two cross-assembly tables."""
    input:
        # Viral contigs are always needed; CheckV and the annotation chain are
        # skipped entirely when annotate is off, matching the databases that
        # Snakefile_phage requires in that mode. The summary tables are still
        # produced, just with the annotation columns left at zero.
        viral = expand(OUTDIR/'phage/{name}/viral/{name}.{caller}.fna',
                       name=sorted(ASSEMBLIES), caller=CALLERS),
        checkv = expand(OUTDIR/'phage/{name}/checkv/{caller}/quality_summary.tsv',
                        name=sorted(ASSEMBLIES), caller=CALLERS) if ANNOTATE else [],
        cds = expand(OUTDIR/'phage/{name}/annotate/{caller}/phynteny/phynteny_per_cds_funcions.tsv',
                     name=sorted(ASSEMBLIES), caller=CALLERS) if ANNOTATE else [],
    output:
        predictions = OUTDIR/'phage/summary/phage_predictions.tsv',
        annotations = OUTDIR/'phage/summary/phage_annotations.tsv'
    params:
        phagedir = OUTDIR/'phage',
        callers = ','.join(CALLERS),
        names = ','.join(sorted(ASSEMBLIES)),
        qerrfile = OUTDIR/'logs/phage/summary.qerr',
        qoutfile = OUTDIR/'logs/phage/summary.qout',
        scratch = 1000,
        mem = 4000,
        time = 60
    log:
        log = OUTDIR/'logs/phage/summary.log'
    threads:
        1
    shell:
        'python ./scripts/phage_summary.py --phage-dir {params.phagedir} '
        '--names {params.names} --callers {params.callers} '
        '--predictions {output.predictions} --annotations {output.annotations} '
        '&> {log.log}'
