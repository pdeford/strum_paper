import glob

def fasta_reader(file_obj):
    """Read a FASTA file and return the sequences in a list."""
    sequences = []
    currSeq = ""
    firstLine = True
    for line in file_obj:
        if line.startswith(">") or line.strip() == "":
            if firstLine:
                firstLine = False
            else:
                sequences.append(currSeq.upper())
                currSeq = ""
        else:
            s = line.strip()
            currSeq += s
        if line == "":
            break
    if currSeq != "":
        sequences.append(currSeq.upper())
    return sequences


assignments = {}

f = open("output/families.txt")
for line in f:
    fields = line.split('\t')
    tf = fields[0]
    family = fields[1]
    if family not in assignments:
        assignments[family] = []
    assignments[family].append(tf)

keepers = []
for key in assignments:
    if len(assignments[key]) >= 10:
        keepers.append(key)

print keepers
for fam in keepers:
    with open('data/{}.fa'.format(fam), 'wb') as g:
        print 'data/{}.fa'.format(fam)
        for tf in assignments[fam]:
            paths = glob.glob('data/{}.ENCFF??????.fa'.format(tf))
            for p in paths:
                accession = p.split('.')[-2]
                with open(p) as f:
                    sequences = fasta_reader(f)
                for i,s in enumerate(sequences[500:1000]):
                    g.write(">{}.{}.{}\n".format(tf, accession, i+500))
                    g.write("{}\n".format(s))
            paths = glob.glob('data/eGFP-{}.ENCFF??????.fa'.format(tf))
            for p in paths:
                accession = p.split('.')[-2]
                with open(p) as f:
                    sequences = fasta_reader(f)
                for i,s in enumerate(sequences[500:1000]):
                    g.write(">{}.{}.{}\n".format(tf, accession, i+500))
                    g.write("{}\n".format(s))
