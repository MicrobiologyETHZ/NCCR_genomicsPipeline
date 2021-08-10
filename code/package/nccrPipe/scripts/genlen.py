# Adapted from Chris Field
# August 01, 2021

import argparse, os, subprocess
from Bio import SeqIO


def blast_genes(faa, db, threads):
    """
    Use diamond to perform sequence alignments of protein coding sequences against a database

    :param faa: fasta file with amino acid sequences
    :param db: DIAMOND database
    :param threads: number of compute threads
    :return: diamond output as a dictionary

    """
    arg = f'diamond blastp --threads {threads} --max-target-seqs 1 --db {db} --query {faa} --outfmt 6 qseqid sseqid slen pident'
    results = subprocess.check_output(arg, shell=True).decode()
    hits = [line.split("\t") for line in results.splitlines()]
    hits = {hit[0]: hit[1:] for hit in hits}
    return hits


def __main__():
    parser = argparse.ArgumentParser(description='Call genes with prodigal, calculate gene length distribution and optionally compare to the nearest BLAST hits')
    parser.add_argument('contigs', metavar='contigs_file', help='File containing contigs in fasta format')
    parser.add_argument('--blast', action='store_true', help='Compare gene lengths to nearest BLAST hits in db')
    parser.add_argument('--db', metavar='ref_db', default='/nfs/cds/Databases/DIAMOND/nr', help='Diamond reference database for gene alignment')
    parser.add_argument('--cutoff', metavar='id_cutoff', default=97.0, help='Identity cutoff for a valid hit')
    parser.add_argument('-t,', '--threads', metavar='threads', default=16, help='Number of compute threads')

    args = parser.parse_args()

    # Determine file names
    contig_dir = os.path.split(args.contigs)[0]
    outdir = os.path.join(contig_dir, "genlen_results")
    prefix = os.path.splitext(os.path.split(args.contigs)[1])[0]
    outfix = os.path.join(outdir, prefix)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    print('Running prodigal')
    # Call genes with prodigal
    prodigal_output = subprocess.check_output(f'prodigal -c -i {args.contigs} -o {outfix}.gff '
                                              f'-a {outfix}.faa -d {outfix}.fna -f gff &> {outfix}.log', shell=True).decode()
    # Read in protein sequences
    proteins = SeqIO.parse(f'{outfix}.faa', 'fasta')
    # Determine gene lengths
    gene_lengths = {protein.id: len(protein) for protein in proteins}

    if args.blast:
        # Blast protein sequences against reference database
        hits = blast_genes(f'{outfix}.faa', args.db, args.threads)
        # Filter hits
        good_hits = {k: v for k, v in hits.items() if float(v[2]) >= args.cutoff}
        # Determine hit lengths
        hit_lengths = {k: int(v[1]) for k, v in good_hits.items()}
        # Calculate ratios
        length_ratios = {k: gene_lengths[k]/hit_lengths[k] for k in hit_lengths.keys()}

    # Output
    with open(f'{outfix}_results.txt', 'w') as fo:
        fo.write('hit\tgene_length\thit_length\tgene_over_hit_length')
        for k in gene_lengths.keys():
            try:
                fo.write(f'{k}\t{gene_lengths[k]}\t{hit_lengths[k]}\t{length_ratios[k]}\n')
            except:
                fo.write(f'{k}\t{gene_lengths[k]}\tNA\tNA\n')

if __name__ == "__main__":
    __main__()

