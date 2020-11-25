import argparse
import subprocess
import shlex
import shutil
import os
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", default='configs/test_config.yaml', type=str,
                        help="Config file")
    parser.add_argument("-a", "--analysis",
                        help="Options: ...", required=True)
    parser.add_argument('-np', '--np', help="Dry run", action='store_true',
                        required=False)
    return parser.parse_args()


def snakemake_cmd(args):
    if args.np:
        cmd = shlex.split(f'snakemake --configfile {args.config} -np {args.analysis} ')
    else:
        rstring = r'"DIR=$(dirname {params.qoutfile}); mkdir -p \"${{DIR}}\"; qsub -S /bin/bash -V -cwd -o {params.qoutfile} -e {params.qerrfile} -pe smp {threads} -l h_vmem={params.mem}M"'
        part1 = shlex.split(f'snakemake --configfile {args.config} --use-conda -k --cluster ')
        part2 = shlex.split(f'{rstring}')
        part3 = shlex.split(f' -p -j 4 --max-jobs-per-second 1 {args.analysis}')
        cmd = part1 + part2 + part3
    return cmd


def main():
    args = parse_args()
    wdPath = Path(__file__).parent.absolute()
    default_config_path = str(wdPath/'configs/test_config.yaml')
    if args.analysis == 'create_config':
        if not os.path.isfile(args.config):
            shutil.copyfile(default_config_path, args.config)
        else:
            print("Config file already exists")
    else:
        cmd = snakemake_cmd(args)
        print(" ".join(cmd))
        subprocess.check_call(cmd, cwd=wdPath)
        print('Done!')


if __name__ == "__main__":
    main()