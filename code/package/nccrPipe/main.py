import argparse
import subprocess
import shlex
import shutil
import os
from pathlib import Path

from .scripts import configure_project

import click

@click.group()
def main():
    pass


#RNAseq
@main.command()
@click.option('--config', '-c', default='configs/rnaseq_config.yaml', help='Configuration File')
@click.option('--method', '-m', default='star', help='star (default) or kallisto')
@click.option('--local',  is_flag=True, help="Run on local machine")
@click.option('--dry',  is_flag=True, help="Show commands without running them")
def rnaseq(config, method, local, dry):
    click.echo("Running Eukaryotic RNASeq Pipeline")
    click.echo(f"Config file: {config}")
    #click.echo("Samples found: ")
    click.echo("Running {}".format('locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile_rnaseq"
    cmd = snakemake_cmd(config, method, smk_file, dry, local)
    click.echo(" ".join(cmd))

# Isolate
@main.command()
@click.option('--config', '-c', default='configs/test_variant_calling_config.yaml', help='Configuration File')
@click.option('--method', '-m', default='call_variants', help='Workflow to run')
@click.option('--local',  is_flag=True, help="Run on local machine")
@click.option('--no-conda',  is_flag=True, help="Do not use conda, under construction")
@click.option('--dry',  is_flag=True, help="Show commands without running them")
def isolate(config, method, local, dry, no_conda):
    click.echo("Running Genomics Pipeline")
    click.echo(f"Config file: {config}")
    #click.echo("Samples found: ")
    click.echo("Running {}".format('locally' if local else ('dry' if dry else 'on cluster')))
    smk_file = "Snakefile"
    cmd = snakemake_cmd(config, method, smk_file, dry, local, no_conda)
    click.echo(" ".join(cmd))


@main.command()
@click.option('--config', '-c', default='configs/basic_config.yaml', help='Unlock working directory if snakemake failed')
def unlock(config):
    cmd = shlex.split(f'snakemake --configfile {config} -j 1 --unlock ')
    wdPath = Path(__file__).parent.absolute()
    subprocess.check_call(cmd, cwd=wdPath)


def snakemake_cmd(config, analysis, smk_file, dry, local, no_conda=False):
    if dry:
        cmd = shlex.split(f'snakemake -s {smk_file} --configfile {config} -np {analysis} ')
    elif local:
        cmd = shlex.split(f'snakemake -s {smk_file} --configfile {config} -j 1 {analysis} ')
    else:
        rstring = r'"DIR=$(dirname {params.qoutfile}); mkdir -p \"${{DIR}}\"; qsub -S /bin/bash -V -cwd -o {params.qoutfile} -e {params.qerrfile} -pe smp {threads} -l h_vmem={params.mem}M"'
        if no_conda:
            part1 = shlex.split(f'snakemake --configfile {config} -s {smk_file} -k --cluster ')
        else:
            part1 = shlex.split(f'snakemake --configfile {config} -s {smk_file} --use-conda -k --cluster ')
        part2 = shlex.split(f'{rstring}')
        part3 = shlex.split(f' -p -j 6 --max-jobs-per-second 1 {analysis}')
        cmd = part1 + part2 + part3
    wdPath = Path(__file__).parent.absolute()
    subprocess.check_call(cmd, cwd=wdPath)
    return cmd

if __name__ == "__main__":
    main()

#
# def parse_args():
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-c", "--config", default='configs/test_config.yaml', type=str,
#                         help="Config file")
#     parser.add_argument("-a", "--analysis",
#                         help="Options: ...", required=True)
#     parser.add_argument('-np', '--np', help="Dry run", action='store_true',
#                         required=False)
#     return parser.parse_args()
#
#
# def snakemake_cmd(args):
#     if args.np:
#         cmd = shlex.split(f'snakemake --configfile {args.config} -np {args.analysis} ')
#     else:
#         rstring = r'"DIR=$(dirname {params.qoutfile}); mkdir -p \"${{DIR}}\"; qsub -S /bin/bash -V -cwd -o {params.qoutfile} -e {params.qerrfile} -pe smp {threads} -l h_vmem={params.mem}M"'
#         part1 = shlex.split(f'snakemake --configfile {args.config} --use-conda -k --cluster ')
#         part2 = shlex.split(f'{rstring}')
#         part3 = shlex.split(f' -p -j 6 --max-jobs-per-second 1 {args.analysis}')
#         cmd = part1 + part2 + part3
#     return cmd
#
#
#
# # def use_modules(args):
# #
# #     cmd_blueprint = f'''
# #     #$ -S /bin/bash
# #     #$ -N {name}
# #     #$ -V
# #     #$ -pe smp {}
# #     #$ -l h_vmem=6G
# #     #$ -e primers.error.log
# #     #$ -o primers.out.log
# #     #$ -t 1-40
# #
# #     {modules}
# #
# #     '''
# #     return cmd
#
#
# def main():
#     args = parse_args()
#     wdPath = Path(__file__).parent.absolute()
#     default_config_path = str(wdPath/'configs/test_config.yaml')
#     if args.analysis == 'create_config':
#         if not os.path.isfile(args.config):
#             #shutil.copyfile(default_config_path, args.config)
#             new_settings = configure_project.configure(default_config_path, args.config)
#         else:
#             print("Config file already exists")
#     else:
#         cmd = snakemake_cmd(args)
#         print(" ".join(cmd))
#         subprocess.check_call(cmd, cwd=wdPath)
#         print('Done!')
#
#
# if __name__ == "__main__":
#     main()