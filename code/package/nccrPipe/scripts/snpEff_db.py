# Add genome entry in SnpEff's snpEff.config
# Add codon table param

# get fasta
# get gbk
from Bio import SeqIO
import sys
'''

ref -> snpEff/data/LL25/genes.gbk

snpEff build -genbank -v LL25


# LL25
LL25.genome: LL25
    LL25.chromosomes : list chromosomes
    LL25.chr.codonTable: Bacterial_and_Plant_Plastid

'''



def write_config_addition(gbk, name, outfile):

    chrs = [seq.id for seq in SeqIO.parse(gbk, "genbank")]
    str1 = f'# {name}\n{name}.genome : {name}\n' \
          f'\t{name}.chromosomes : {", ".join(chrs)}\n'
    str2 = ''
    for c in chrs:
        str2+=f'\t{name}.{c}.codonTable : Bacterial_and_Plant_Plastid\n'
    with open(outfile, 'w') as fo:
        fo.write(str1+str2)



def rename_vcf(vcf_file_in, vcf_file_out, genome='LL25'):
    with open(vcf_file_in, 'r') as fh:
        with open (vcf_file_out, 'w') as fo:
            for line in fh.readlines():
                w = line.split()
                N = w[0].split('_')[1]
                newN = f'{genome}_{N}'
                fo.write("\t".join([newN] + w[1:]) + "\n")


# write_config_addition("/science/ansintsova/esbl_strains/leo_project2/pair1/annotation/LL25/LL25.gbk",
#                       'LL25', "/nfs/home/ansintsova/esbl_strains/LL25.snpEff.txt")


def add_to_snpEff(env = "/nfs/home/ansintsova/./miniconda3/pkgs/snpeff-5.0-0/share/snpeff-5.0-0/"): # todo not the right env (right only for config)
    #Create text file
    # append to config

    # Copy gbk to env/data/GENOME/genes.gbk
    # run snpEff build -genbank -v GENOME
    return None

if __name__ == "__main__":
    try:
        if sys.argv[1] == 'help':
            print('Option: rename (rename VCF with ref conitgs), add (create a file to be appened to snpEff config')
        elif sys.argv[1] == 'rename':

            assert len(sys.argv) > 4
            vcfIn = sys.argv[2]
            vcfOut = sys.argv[3]
            genome = sys.argv[4]
            rename_vcf(vcfIn, vcfOut, genome)
        elif sys.argv[1] == 'add':
            assert len(sys.argv) > 4
            gbk = sys.argv[2]
            name = sys.argv[3]
            outfile = sys.argv[4]
            write_config_addition(gbk, name, outfile)
        print(sys.argv[1])
    except IndexError:
        print("Provide an argument, help for options")
    #rename_vcf('/nfs/home/ansintsova/esbl_strains/LL28.test.vcf',
    # '/nfs/home/ansintsova/esbl_strains/LL28.rename.vcf', 'LL25')

