import argparse
import importlib.util
import subprocess
import shlex
import shutil
import os
from pathlib import Path

from .scripts import configure_project

import click
from .scripts import fastq_dir_to_samplesheet as fds
import yaml


@click.group()
def main():
    pass


@main.command(help="Generate samplesheet from a directory of FastQ files.")
@click.option('--configfile', '-c', default='', help='Configuration File')
@click.option("-i", "--fastq_dir", help="Folder containing raw FastQ files.")
@click.option("-o", "--sample_file", help="Output samplesheet file.")
@click.option("-r1", "--read1_extension", type=str,  default="_R1.fq.gz",
              help="File extension for read 1.")
@click.option("-r2", "--read2_extension", type=str, default="_R2.fq.gz",
              help="File extension for read 2.")
@click.option("-sn", "--sanitise_name", is_flag=True,
              help="Whether to further sanitise FastQ file name to get sample id. Used in conjunction with "
                   "--sanitise_name_delimiter and --sanitise_name_index.")
@click.option("-sd", "--sanitise_name_delimiter", type=str, default="_",
              help="Delimiter to use to sanitise sample name.", )
@click.option("-si", "--sanitise_name_index", type=int, default=1,
              help="After splitting FastQ file name by --sanitise_name_delimiter "
                   "all elements before this index (1-based) will be joined to create final sample name.",)
def samples(configfile, fastq_dir, sample_file, read2_extension, read1_extension, sanitise_name,
            sanitise_name_delimiter, sanitise_name_index):
    click.echo("Generating samples file from fastq directory...")
    if configfile:
        click.echo(f"Config file: {configfile}")
        with open(configfile) as file:
            config = yaml.load(file, Loader=yaml.FullLoader)
        # Paths in configs are relative to workflow/ (where snakemake runs).
        # Resolve them to absolute so glob works regardless of the user's CWD.
        workflow_dir = Path(__file__).parent
        def resolve_path(p):
            p = Path(p)
            return p if p.is_absolute() else (workflow_dir / p).resolve()
        abs_fastq_dir = resolve_path(config['dataDir'])
        abs_samplesheet = resolve_path(config['samples'])
        fds.fastq_dir_to_samplesheet(
            fastq_dir=str(abs_fastq_dir),
            samplesheet_file=str(abs_samplesheet),
            read1_extension=config['fq_fwd'],
            read2_extension=config['fq_rvr'],
            sanitise_name=config['sanitise_name'],
            sanitise_name_delimiter=config['name_delimiter'],
            sanitise_name_index=config['name_index'],
        )
    else:
        fds.fastq_dir_to_samplesheet(
            fastq_dir=fastq_dir,
            samplesheet_file=sample_file,
            read1_extension=read1_extension,
            read2_extension=read2_extension,
            sanitise_name=sanitise_name,
            sanitise_name_delimiter=sanitise_name_delimiter,
            sanitise_name_index=sanitise_name_index,
        )


PARTITION_OPTION = click.option('--partition', '-p', default='institute', show_default=True,
                               help="SLURM partition to use for cluster jobs")

# Preprocess
@main.command()
@click.option('--config', '-c', default='configs/test_variant_calling_config.yaml', help='Configuration File')
@click.option('--local',  is_flag=True, help="Run on local machine")
@click.option('--no-conda',  is_flag=True, help="Do not use conda, under construction")
@click.option('--dry',  is_flag=True, help="Show commands without running them")
@PARTITION_OPTION
def clean(config, local, dry, no_conda, partition):
    click.echo("Running Preprocessing Pipeline")
    click.echo(f"Config file: {config}")
    click.echo("Running {}".format(
        'locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile"
    cmd = snakemake_cmd(config, 'preprocess', smk_file, dry, local, no_conda, partition)
    click.echo(" ".join(cmd))


# Assembly
@main.command()
@click.option('--config', '-c', default='configs/test_variant_calling_config.yaml', help='Configuration File')
@click.option('--local',  is_flag=True, help="Run on local machine")
@click.option('--no-conda',  is_flag=True, help="Do not use conda, under construction")
@click.option('--dry',  is_flag=True, help="Show commands without running them")
@PARTITION_OPTION
def assemble(config, local, dry, no_conda, partition):
    click.echo("Running Assembly Pipeline")
    click.echo(f"Config file: {config}")
    click.echo("Running {}".format(
        'locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile"
    cmd = snakemake_cmd(config, 'assemble', smk_file, dry, local, no_conda, partition)
    click.echo(" ".join(cmd))


# Align
@main.command()
@click.option('--config', '-c', default='configs/test_variant_calling_config.yaml', help='Configuration File')
@click.option('--local',  is_flag=True, help="Run on local machine")
@click.option('--no-conda',  is_flag=True, help="Do not use conda, under construction")
@click.option('--dry',  is_flag=True, help="Show commands without running them")
@PARTITION_OPTION
def align(config, local, dry, no_conda, partition):
    click.echo("Running Align Command")
    click.echo(f"Config file: {config}")
    click.echo("Running {}".format(
        'locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile"
    cmd = snakemake_cmd(config, 'align', smk_file, dry, local, no_conda, partition)
    click.echo(" ".join(cmd))


# Call
@main.command()
@click.option('--config', '-c', default='configs/test_variant_calling_config.yaml', help='Configuration File')
@click.option('--local',  is_flag=True, help="Run on local machine")
@click.option('--no-conda',  is_flag=True, help="Do not use conda, under construction")
@click.option('--dry',  is_flag=True, help="Show commands without running them")
@PARTITION_OPTION
def call(config, local, dry, no_conda, partition):
    click.echo("Running Breseq Variant Calling Pipeline")
    click.echo(f"Config file: {config}")
    click.echo("Running {}".format(
        'locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile"
    cmd = snakemake_cmd(config, 'breseq', smk_file, dry, local, no_conda, partition)
    click.echo(" ".join(cmd))


# Call
@main.command()
@click.option('--config', '-c', default='configs/test_variant_calling_config.yaml', help='Configuration File')
@click.option('--local',  is_flag=True, help="Run on local machine")
@click.option('--no-conda',  is_flag=True, help="Do not use conda, under construction")
@click.option('--dry',  is_flag=True, help="Show commands without running them")
@PARTITION_OPTION
def funcall(config, local, dry, no_conda, partition):
    click.echo("Running Eukaryotic Variant Calling Pipeline")
    click.echo(f"Config file: {config}")
    click.echo("Running {}".format(
        'locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile"
    cmd = snakemake_cmd(config, 'varcall', smk_file, dry, local, no_conda, partition)
    click.echo(" ".join(cmd))


# InStrain
@main.command()
@click.option('--config', '-c', default='configs/instrain_config.yaml', help='Configuration File')
@click.option('--compare/--no-compare', default=False,
              help="Also run inStrain compare across samples (target: instrain_compare)")
@click.option('--local',  is_flag=True, help="Run on local machine")
@click.option('--no-conda',  is_flag=True, help="Do not use conda, under construction")
@click.option('--dry',  is_flag=True, help="Show commands without running them")
@PARTITION_OPTION
def instrain(config, compare, local, dry, no_conda, partition):
    click.echo("Running InStrain Microdiversity Pipeline")
    click.echo(f"Config file: {config}")
    click.echo("Running {}".format(
        'locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile"
    target = 'instrain_compare' if compare else 'instrain'
    cmd = snakemake_cmd(config, target, smk_file, dry, local, no_conda, partition)
    click.echo(" ".join(cmd))


# ISMapper
@main.command()
@click.option('--config', '-c', default='configs/ismap_config.yaml', help='Configuration File')
@click.option('--local',  is_flag=True, help="Run on local machine")
@click.option('--no-conda',  is_flag=True, help="Do not use conda, under construction")
@click.option('--dry',  is_flag=True, help="Show commands without running them")
@PARTITION_OPTION
def ismap(config, local, dry, no_conda, partition):
    """Locate insertion sequence (IS) sites with ISMapper."""
    click.echo("Running ISMapper Insertion Sequence Pipeline")
    click.echo(f"Config file: {config}")
    click.echo("Running {}".format(
        'locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile"
    cmd = snakemake_cmd(config, 'ismap', smk_file, dry, local, no_conda, partition)
    click.echo(" ".join(cmd))


# Annotate
@main.command()
@click.option('--config', '-c', default='configs/test_variant_calling_config.yaml', help='Configuration File')
@click.option('--local',  is_flag=True, help="Run on local machine")
@click.option('--no-conda',  is_flag=True, help="Do not use conda, under construction")
@click.option('--dry',  is_flag=True, help="Show commands without running them")
@PARTITION_OPTION
def annotate(config, local, dry, no_conda, partition):
    click.echo("Running Assembly Pipeline")
    click.echo(f"Config file: {config}")
    click.echo("Running {}".format(
        'locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile"
    cmd = snakemake_cmd(config, 'annotate', smk_file, dry, local, no_conda, partition)
    click.echo(" ".join(cmd))


# New and exploratory
@main.command()
@click.option('--config', '-c', default='configs/test_variant_calling_config.yaml', help='Configuration File')
@click.option('--local',  is_flag=True, help="Run on local machine")
@click.option('--no-conda',  is_flag=True, help="Do not use conda, under construction")
@click.option('--dry',  is_flag=True, help="Show commands without running them")
@PARTITION_OPTION
def gapseq(config, local, dry, no_conda, partition):
    click.echo("Running Gapseq find Pipeline")
    click.echo(f"Config file: {config}")
    # click.echo("Samples found: ")
    click.echo("Running {}".format(
        'locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile"
    cmd = snakemake_cmd(config, 'run_gapseq', smk_file, dry, local, no_conda, partition)
    click.echo(" ".join(cmd))


# Isolate
@main.command()
@click.option('--config', '-c', default='configs/test_variant_calling_config.yaml', help='Configuration File')
@click.option('--method', '-m', default='call_variants', help='Workflow to run, '
                                                              'options: [call_variants, assemble, assemble_only]')
@click.option('--local',  is_flag=True, help="Run on local machine")
@click.option('--no-conda',  is_flag=True, help="Do not use conda, under construction")
@click.option('--dry',  is_flag=True, help="Show commands without running them")
@PARTITION_OPTION
def isolate(config, method, local, dry, no_conda, partition):
    click.echo("Running Genomics Pipeline - UNDER CONSTRUCTION - use assemble")
    click.echo(f"Config file: {config}")
    # click.echo("Samples found: ")
    click.echo("Running {}".format(
        'locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile"
    cmd = snakemake_cmd(config, method, smk_file, dry, local, no_conda, partition)
    click.echo(" ".join(cmd))


@main.command()
@click.option('--config', '-c', default='configs/basic_config.yaml', help='Unlock working directory if snakemake failed')
def unlock(config):
    cmd = shlex.split(f'snakemake --configfile {config} -j 1 --unlock ')
    wdPath = Path(__file__).parent.absolute()
    subprocess.check_call(cmd, cwd=wdPath)


# TODO: under construction

@main.command()
@click.option('-c', '--config', required=True, help='Configuration file')
@click.option('-m', '--method', default='metaflye', help='Workflow to run [metaflye]')
@click.option('--local', is_flag=True, help='Run locally')
@click.option('--dry', is_flag=True, help='Dry run')
@PARTITION_OPTION
def metagenome(config, method, local, dry, partition):
    """Run metagenomic assembly comparison"""
    click.echo("Running Metagenomic Workflow")
    click.echo(f"Running {method}")
    smk_file = Path(__file__).parent / "Snakefile_metagenome"
    cmd = snakemake_cmd(config, method, smk_file, dry, local, partition=partition)
    click.echo(" ".join(cmd))


@main.command()
@click.option('-c', '--config', required=True, help='Configuration file')
@click.option('--predict-only', is_flag=True,
              help="Stop after prediction; skip CheckV and the annotation chain")
@click.option('--summary/--no-summary', default=True,
              help="Also build the merged cross-assembly summary tables")
@click.option('--local', is_flag=True, help='Run locally')
@click.option('--no-conda', is_flag=True, help="Do not use conda, under construction")
@click.option('--dry', is_flag=True, help='Dry run')
@click.option('--cores', '-j', type=int, default=None,
              help="Cores for a local run / concurrent jobs on the cluster")
@click.option('--set', 'set_config', multiple=True, metavar='KEY=VALUE',
              help="Override a config value, e.g. --set annotate=false. Repeatable.")
@click.option('--latency-wait', type=int, default=60, show_default=True,
              help="Seconds to wait for output files to appear on a shared "
                   "filesystem after a cluster job finishes")
@PARTITION_OPTION
def phage(config, predict_only, summary, local, dry, no_conda, cores,
          set_config, latency_wait, partition):
    """Predict and annotate phages in assemblies.

    Prediction with geNomad and Cenote-Taker 3, annotation with
    pharokka -> phold -> phynteny_transformer, quality with CheckV.

    \b
    Examples:
      # Dry run: show the DAG without running anything
      nccrPipe phage -c configs/test_phage_config.yaml --dry --predict-only

      # Smallest real run: one caller, prediction only, 8 cores
      nccrPipe phage -c my_config.yaml --local -j 8 --predict-only \\
          --set callers=genomad

      # Full prediction + annotation + summary tables on the cluster
      nccrPipe phage -c my_config.yaml
    """
    if predict_only:
        target = 'find_phage'
    elif summary:
        target = 'phage_summary'
    else:
        target = 'annotate_phage'
    click.echo("Running phage prediction/annotation workflow")
    click.echo(f"Config file: {config}")
    click.echo(f"Target: {target}")
    click.echo("Running {}".format(
        'locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = Path(__file__).parent / "Snakefile_phage"
    cmd = snakemake_cmd(config, target, smk_file, dry, local, no_conda,
                        partition, cores=cores, set_config=set_config,
                        latency_wait=latency_wait)
    click.echo(" ".join(cmd))


def snakemake_cmd(config, analysis, smk_file, dry, local, no_conda=False,
                  partition='institute', cores=None, set_config=(),
                  latency_wait=60):
    config_path = Path(config)
    if not config_path.is_absolute():
        resolved = config_path.resolve()
        if not resolved.exists():
            repo_root = Path(__file__).parent.parent
            resolved = repo_root / config_path
        config_path = resolved
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config}")
    config = str(config_path)
    # `--config k=v` overrides, appended last so they win over the config file.
    overrides = ['--config', *set_config] if set_config else []
    if dry:
        cmd = shlex.split(
            f'snakemake -s {smk_file} --configfile {config} -np {analysis} ')
    elif local:
        # Honour --no-conda here as well as on the cluster; without --use-conda
        # the `conda:` directives are ignored entirely and every tool has to be
        # on the PATH of the active environment.
        conda_arg = '' if no_conda else '--use-conda '
        cmd = shlex.split(
            f'snakemake -s {smk_file} --configfile {config} {conda_arg}'
            f'-j {cores or 1} {analysis} ')
    else:
        check_cluster_executor()
        conda_arg = '' if no_conda else '--use-conda '
        # The submit command must reach snakemake as a SINGLE argv element, so it
        # is appended directly rather than round-tripped through shlex.split.
        # The {params.*} / {threads} placeholders are expanded by snakemake per
        # job, not by us, so they must survive untouched.
        cmd = shlex.split(
            f'snakemake --configfile {config} -s {smk_file} {conda_arg}-k '
            f'--executor cluster-generic --cluster-generic-submit-cmd ')
        cmd.append(slurm_submit_cmd(partition))
        cmd += shlex.split(
            f'-p -j {cores or 6} --max-jobs-per-second 1 '
            f'--latency-wait {latency_wait} {analysis}')
    cmd += overrides
    wdPath = Path(__file__).parent.absolute()
    subprocess.check_call(cmd, cwd=wdPath)
    return cmd


def slurm_submit_cmd(partition):
    """Build the sbatch command that cluster-generic runs for each job.

    Resources come from each rule's `params:` (mem, time) and `threads`. The
    braces are snakemake placeholders expanded per job, so this string is
    deliberately not formatted here beyond the partition.
    """
    return (
        'DIR=$(dirname {params.qoutfile}); mkdir -p "$DIR"; '
        'sbatch -t {params.time} --mem-per-cpu={params.mem} -n {threads} '
        '-o {params.qoutfile} -e {params.qerrfile} '
        f'--partition {partition}'
    )


def check_cluster_executor():
    """Fail early and actionably if the cluster executor plugin is missing.

    Snakemake 8 deprecated and 9 removed `--cluster`; submission now needs an
    executor plugin. Without this check the failure is snakemake's
    `unrecognized arguments: --cluster-generic-submit-cmd`, which gives no hint
    as to what is actually wrong.
    """
    if importlib.util.find_spec('snakemake_executor_plugin_cluster_generic') is None:
        raise click.ClickException(
            "Cluster submission needs the cluster-generic executor plugin, which "
            "is not installed.\n\n"
            "    conda install -c conda-forge -c bioconda "
            "snakemake-executor-plugin-cluster-generic\n\n"
            "(Snakemake 9 removed the old --cluster option.) "
            "Alternatively run with --local, or preview with --dry."
        )


if __name__ == "__main__":
    main()

