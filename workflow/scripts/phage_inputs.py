"""
Input resolution for the phage detection / annotation workflow.

The phage workflow is usually pointed at assemblies this pipeline did *not*
produce (public genomes, collaborator data), so it cannot be driven by
samples.csv alone. `resolve_assemblies` collects assemblies from any mix of
three config sources and returns a ``{name: absolute_path}`` mapping. Every
rule in rules/phage.smk is keyed on that short name rather than on the input
path, which keeps outputs under OUTDIR and off read-only input directories.

Config sources (all optional, merged in this order so later ones win):

    assembly_dir: /path/to/assemblies    # + optional `pattern:` (str or list)
    assemblies:                          # list of paths ...
      - /path/to/asm1.fna
    assemblies:                          # ... or an explicit name -> path map
      my_genome: /path/to/asm1.fna
    samples: samples.csv                 # assemblies built by this pipeline

This module is plain Python with no Snakemake imports so it can be unit tested
directly; see tests/unit/test_phage_inputs.py.
"""
import sys
from pathlib import Path

# Suffixes stripped when deriving a name from a filename. Compression suffixes
# are stripped first, so `GCA_000001_genomic.fna.gz` -> `GCA_000001_genomic`.
COMPRESSION_SUFFIXES = {'.gz', '.bz2', '.xz', '.zst'}
FASTA_SUFFIXES = {'.fa', '.fasta', '.fna', '.fas', '.ffn', '.contigs'}

# Filename patterns globbed under `assembly_dir` when `pattern` is not given.
DEFAULT_PATTERNS = ['*.fa', '*.fasta', '*.fna', '*.fas']

# Databases the workflow needs, with the command that installs each one. Used to
# produce an actionable error at parse time instead of a tool failure an hour in.
DATABASE_INSTALL = {
    'genomad': 'genomad download-database <dir>',
    'cenotetaker': ('get_ct3_dbs -o <dir> --hmm T --hallmark_tax T '
                    '--refseq_tax T --mmseqs_cdd T --domain_list T'),
    'checkv': 'checkv download_database <dir>',
    'pharokka': 'install_databases.py -o <dir>',
    'phold': 'phold install -d <dir>',
    'phynteny': 'install_models -o <dir>',
}


# Values that mean "off" when a flag arrives as a string. Snakemake's
# `--config key=false` does NOT yield a bool — it leaves the string "false",
# which is truthy — so every boolean config flag must go through as_bool.
FALSEY = {'false', 'no', 'off', '0', 'none', ''}


def as_bool(value, default=True):
    """Interpret a config value as a boolean, tolerating strings from --config."""
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() not in FALSEY
    return bool(value)


def resolve_path(p, basedir):
    """Resolve a config path to absolute, relative to the workflow directory.

    Mirrors the `_resolve` helper in workflow/Snakefile so that paths in phage
    configs behave the same as paths in the main pipeline's configs.
    """
    p = Path(p)
    return p if p.is_absolute() else (Path(basedir) / p).resolve()


def assembly_name(path):
    """Derive a short name from an assembly filename.

    Strips one compression suffix and one FASTA suffix, so both
    `scaffolds.fasta` and `GCA_000001_genomic.fna.gz` yield sensible names.
    """
    p = Path(path)
    if p.suffix.lower() in COMPRESSION_SUFFIXES:
        p = p.with_suffix('')
    if p.suffix.lower() in FASTA_SUFFIXES:
        p = p.with_suffix('')
    return p.name


def is_compressed(path):
    """True if the assembly needs decompressing before the tools can read it."""
    return Path(path).suffix.lower() in COMPRESSION_SUFFIXES


def _glob_assembly_dir(assembly_dir, pattern):
    """Collect assembly files under `assembly_dir` matching `pattern`.

    `pattern` may be a single glob or a list of globs. Each pattern also matches
    its compressed form, since public assemblies almost always arrive gzipped.
    """
    if pattern is None:
        patterns = list(DEFAULT_PATTERNS)
    elif isinstance(pattern, str):
        patterns = [pattern]
    else:
        patterns = list(pattern)

    # Match `*.fna` and `*.fna.gz` alike without requiring the user to spell out
    # both. If a caller explicitly asks for a compressed pattern, leave it be.
    expanded = []
    for pat in patterns:
        expanded.append(pat)
        if Path(pat).suffix.lower() not in COMPRESSION_SUFFIXES:
            expanded.extend(f'{pat}{ext}' for ext in sorted(COMPRESSION_SUFFIXES))

    found = set()
    for pat in expanded:
        found.update(p for p in assembly_dir.glob(pat) if p.is_file())
    return sorted(found)


def _assign_names(paths):
    """Map assembly paths to unique short names.

    Uses the filename stem; where that collides, falls back to
    `<parent_dir>_<stem>`. Raises if names still collide after the fallback,
    rather than silently letting one assembly's results overwrite another's.
    """
    by_stem = {}
    for path in paths:
        by_stem.setdefault(assembly_name(path), []).append(path)

    resolved = {}
    for stem, group in by_stem.items():
        if len(group) == 1:
            resolved[stem] = group[0]
            continue
        # Collision (e.g. several pipeline runs all producing `scaffolds.fasta`)
        # — disambiguate with the containing directory.
        for path in group:
            resolved[f'{path.parent.name}_{stem}'] = path

    if len(resolved) != len(paths):
        counts = {}
        for path in paths:
            counts.setdefault(assembly_name(path), []).append(str(path))
        clashes = {stem: p for stem, p in counts.items() if len(p) > 1}
        raise ValueError(
            "Could not derive unique names for these assemblies even after "
            "falling back to <parent_dir>_<filename>:\n"
            + "\n".join(f'  {stem}: {", ".join(paths_)}'
                        for stem, paths_ in clashes.items())
            + "\nGive them explicit names using the mapping form of "
              "`assemblies:` in the config, e.g.\n"
              "  assemblies:\n    my_name: /path/to/scaffolds.fasta"
        )
    return resolved


def _from_samples(config, basedir):
    """Assemblies produced by this pipeline, for the samples in samples.csv.

    Mirrors the assembly output paths built in workflow/Snakefile:96-101.
    """
    import pandas as pd

    samples_file = resolve_path(config['samples'], basedir)
    samples = pd.read_csv(samples_file, comment='#')['sample'].unique()
    outdir = resolve_path(config.get('outDir', config.get('output_dir', 'output')),
                          basedir)
    assembler = config.get('assembler', 'spades')

    if assembler == 'spades':
        return {s: outdir / f'assembly/{s}/{s}.scaffolds.min200.fasta'
                for s in samples}
    if assembler == 'unicycler':
        mode = config.get('unimode', '')
        return {s: outdir / f'unicycler/{mode}/{s}/assembly.fasta' for s in samples}
    raise ValueError(
        f"Unknown assembler '{assembler}'; expected 'spades' or 'unicycler'."
    )


def resolve_assemblies(config, basedir):
    """Build the ``{name: absolute_path}`` mapping driving the phage workflow.

    `basedir` is the workflow directory (``workflow.basedir`` in a Snakefile);
    relative config paths are resolved against it.
    """
    assemblies = {}

    assembly_dir = config.get('assembly_dir', '')
    if assembly_dir:
        assembly_dir = resolve_path(assembly_dir, basedir)
        if not assembly_dir.is_dir():
            raise ValueError(f"assembly_dir does not exist: {assembly_dir}")
        found = _glob_assembly_dir(assembly_dir, config.get('pattern'))
        if not found:
            raise ValueError(
                f"No assemblies found in {assembly_dir} matching "
                f"{config.get('pattern') or DEFAULT_PATTERNS}."
            )
        assemblies.update(_assign_names(found))

    explicit = config.get('assemblies', [])
    if isinstance(explicit, dict):
        # Explicit name -> path mapping always wins; no name derivation needed.
        assemblies.update({name: resolve_path(p, basedir)
                           for name, p in explicit.items()})
    elif explicit:
        assemblies.update(
            _assign_names([resolve_path(p, basedir) for p in explicit]))

    if config.get('samples'):
        assemblies.update(_from_samples(config, basedir))

    if not assemblies:
        raise ValueError(
            "No assemblies to process. Provide at least one of `assembly_dir:`, "
            "`assemblies:` or `samples:` in the phage config."
        )
    return assemblies


def _is_placeholder(path):
    """True if `path` is one of the repo's empty test-placeholder directories."""
    if not path.is_dir():
        return False
    entries = list(path.iterdir())
    return len(entries) == 1 and entries[0].name == '.placeholder'


def validate_databases(config, required, basedir):
    """Check that each required tool's database path is set and non-empty.

    Databases are installed out of band (see workflow/scripts/setup_phage_dbs.sh);
    this only verifies the config points somewhere real, so a 20 GB-download
    mistake surfaces at parse time rather than mid-run.
    """
    databases = config.get('databases', {})
    problems = []
    for tool in sorted(required):
        path = databases.get(tool, '')
        if not path:
            problems.append(
                f"  {tool}: no path set (databases: {tool}: <dir>)\n"
                f"      install with: {DATABASE_INSTALL[tool]}")
            continue
        resolved = resolve_path(path, basedir)
        if not resolved.exists():
            problems.append(
                f"  {tool}: path does not exist: {resolved}\n"
                f"      install with: {DATABASE_INSTALL[tool]}")
        elif _is_placeholder(resolved):
            # configs/test_phage_config.yaml points at empty placeholder dirs so
            # the dry-run tests need no real databases. Those paths get copied
            # into real configs, where they pass an existence check and the
            # failure only surfaces later as an opaque error from inside the
            # tool. Warn rather than raise: the test config legitimately uses
            # these, so erroring here would break every dry run.
            print(
                f"WARNING: {tool} database {resolved} is a repo test "
                f"placeholder, not a real database. Jobs using it will fail.\n"
                f"         install a real one with: {DATABASE_INSTALL[tool]}",
                file=sys.stderr)

    if problems:
        raise ValueError(
            "Missing databases required by the phage workflow:\n"
            + "\n".join(problems)
            + "\n\nAll six can be installed with "
              "workflow/scripts/setup_phage_dbs.sh <parent_dir>."
        )
    return {tool: str(resolve_path(databases[tool], basedir)) for tool in required}
