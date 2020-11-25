file = '/Users/ansintsova/git_repos/TNSEQ_EAN/data/expectedQuantEd.tsv'

fileO = "/Users/ansintsova/git_repos/TNSEQ_EAN/data/expectedQuantEd2.tsv"
with open(file, 'r') as fh:
    with open(fileO, 'w') as fo:
        for line in fh.readlines():
            rep = int(line.split()[1])
            fo.write(line*rep)


