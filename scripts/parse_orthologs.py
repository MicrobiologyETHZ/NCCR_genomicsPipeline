from pathlib import Path
import pandas as pd
import sys

def parseOrthFile(file):
    with open(file, 'r') as fh:
        #genome_to_genes = {}
        genes_to_genome = {}
        geneID = Path(file).stem
        for line in fh.readlines():
            if line.startswith(">"):
                genome = line.strip(">").split("|")[0]
                gene = line.split("|")[1].split("-")[0]
                #genome_to_genes[genome] =  {**genome_to_genes.get(genome, {}), **{geneID:gene}}
                genes_to_genome[geneID] = {**genes_to_genome.get(geneID, {}), **{genome:gene}}
        return  genes_to_genome

def parseAllFiles(files_list, outFile):
    genes = {}
    for file in files_list:
        to_genome = parseOrthFile(file)
        genes = {**genes, **to_genome}
    df = pd.DataFrame(genes).T
    df.to_csv(outFile)
    return df

#################
def parseEmapperAnnotation(files, outFile):
    genes_to_genome = {}
    for file in files:
        with open(file, 'r') as fh:
            for line in fh.readlines():
                if not line.startswith('#'):
                    geneID = line.split('\t')[1]
                    genome = line.split('\t')[0].split("_")[0]
                    gene = line.split()[0].split("_")[1]
                    genes_to_genome[geneID] = {**genes_to_genome.get(geneID, {}), **{genome: gene}}
    df = pd.DataFrame(genes_to_genome).T
    df.to_csv(outFile)
    return df



if __name__ == "__main__":
    file_out = sys.argv[1]
    files = [sys.argv[i] for i in range(2, len(sys.argv))]
    #df = parseAllFiles(file_list, file_out)
    df = parseEmapperAnnotation(files, file_out)
    print(df.head())
