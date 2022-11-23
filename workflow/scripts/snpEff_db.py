# Add genome entry in SnpEff's snpEff.config
# Add codon table param

# get fasta
# get gbk
from Bio import SeqIO
import sys
from pathlib import Path
import shutil
import subprocess
import shlex
import argparse
import os
'''

ref -> snpEff/data/LL25/genes.gbk

snpEff build -genbank -v LL25


# LL25
LL25.genome: LL25
    LL25.chromosomes : list chromosomes
    LL25.chr.codonTable: Bacterial_and_Plant_Plastid

'''

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("gbk")
    parser.add_argument("genome_name")
    parser.add_argument("-nenv",  required=False)
    parser.add_argument("-fa", required=False)
    parser.add_argument('-env',
                        required=False)
    return parser.parse_args()

def config_addition(gbk, name):

    chrs = [seq.id for seq in SeqIO.parse(gbk, "genbank")]
    str1 = f'# {name}\n{name}.genome : {name}\n' \
          f'\t{name}.chromosomes : {", ".join(chrs)}\n'
    str2 = ''
    for c in chrs:
        str2 += f'\t{name}.{c}.codonTable : Bacterial_and_Plant_Plastid\n'
    return str1+str2


def add_to_snpEff(gbk,  genome_name, fasta='', nenv='',
                  env="/nfs/home/ansintsova/./miniconda3/pkgs/snpeff-5.1-hdfd78af_2/share/snpeff-5.1-2/"):
    # Create genome info and append to snpeff.config
    str_to_append = config_addition(gbk, genome_name)
    config_path = Path(env)/"snpEff.config"
    with open(config_path, "a") as fh:
        fh.write(str_to_append)

    # Copy files into appropriate locations and add symlink if needed
    os.makedirs(os.path.dirname(env + f'/data/{genome_name}/genes.gbk'), exist_ok=True)
    shutil.copyfile(gbk, env + f'/data/{genome_name}/genes.gbk')
    print(os.path.isfile(env + f'/data/{genome_name}/genes.gbk'))
    if fasta:
        shutil.copyfile(fasta, env + f'/data/{genome_name}/sequences.fa')
    if nenv:
        subprocess.call(shlex.split(f"ln -s {env}/data {nenv}"))
        os.remove(f'{nenv}/snpEff.config')
        subprocess.call(shlex.split(f"ln -s {env}/snpEff.config {nenv}"))

"""
salmonella_enterica_SL1344_FQ312003.1.gbk
salmonella_enterica_SL1344_FQ312003.1.fasta

"""


if __name__ == "__main__":
    args = parse_args()
    fa = args.fa if args.fa else ''
    nenv = args.nenv if args.nenv else ''
    env = args.env if args.env else '/nfs/home/ansintsova/./miniconda3/pkgs/snpeff-5.1-hdfd78af_2/share/snpeff-5.1-2/'
    add_to_snpEff(args.gbk, args.genome_name, fa, nenv, env)





    """

    # run snpEff build -genbank -v GENOME
    return None



# def rename_vcf(vcf_file_in, vcf_file_out, genome='LL25'):
#     with open(vcf_file_in, 'r') as fh:
#         with open (vcf_file_out, 'w') as fo:
#             for line in fh.readlines():
#                 w = line.split()
#                 N = w[0].split('_')[1]
#                 newN = f'{genome}_{N}'
#                 fo.write("\t".join([newN] + w[1:]) + "\n")


# write_config_addition("/science/ansintsova/esbl_strains/leo_project2/pair1/annotation/LL25/LL25.gbk",
#                       'LL25', "/nfs/home/ansintsova/esbl_strains/LL25.snpEff.txt")


 """