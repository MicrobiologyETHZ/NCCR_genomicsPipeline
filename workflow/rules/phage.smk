"""
Phage / provirus detection with geNomad.

Runs `genomad end-to-end` on each assembly FASTA listed in the config. Unlike the
main pipeline (which is driven by samples.csv), this step takes an explicit list of
assembly FASTA paths via `sample:` — see configs/phage_config.yaml. geNomad outputs
(virus/plasmid summaries, proviruses, etc.) are written next to each input assembly,
with a `<assembly>.genomad.done` marker on success.

Requires a geNomad database (download once with `genomad download-database <dir>`).

Config (see configs/phage_config.yaml):

    sample:                       # list of assembly FASTA paths
      - /path/to/asm1/scaffolds.fasta
    genomad:
      db: /path/to/genomad_db     # geNomad database directory
      disable_nn: true            # --disable-nn-classification (faster; skips NN classifier)
"""
from pathlib import Path

_genomad_cfg = config.get('genomad', {})
# Accept the legacy top-level `genomad_db` key as a fallback for the db path.
GENOMAD_DB = _genomad_cfg.get('db', config.get('genomad_db', ''))
GENOMAD_DISABLE_NN = _genomad_cfg.get('disable_nn', True)


rule run_genomad:
    input:
        assembly = '{assembly}'
    output:
        marker = touch('{assembly}.genomad.done')
    params:
        outdir = lambda wildcards: str(Path(wildcards.assembly).parent),
        db = GENOMAD_DB,
        nn_arg = '--disable-nn-classification' if GENOMAD_DISABLE_NN else '',
        scratch = 1000,
        mem = 4000,
        time = 800,
        qerrfile = lambda wildcards: str(Path(wildcards.assembly).parent.stem) + '.genomad.qerr',
        qoutfile = lambda wildcards: str(Path(wildcards.assembly).parent.stem) + '.genomad.qout'
    conda:
        'phage'
    log:
        log = '{assembly}.genomad.log'
    threads:
        16
    shell:
        'genomad end-to-end --threads {threads} {params.nn_arg} '
        '{input.assembly} {params.outdir} {params.db} &> {log.log}'
