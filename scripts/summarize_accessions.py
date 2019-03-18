#!/usr/bin/env python

outfiles = [
    "output/coefficents.txt",
    "output/alignments.txt",
    "output/chip_auc.txt",
    "output/correlations.txt",
    "output/position_comp.txt",
    "output/specificities.txt",
]

labels = [
    'coefficients',
    'alignments',
    'classification',
    'correlation',
    'positions',
    'specificities',
]

gene_list = []
for outfile in outfiles:
    with open(outfile) as f:
        genes = []
        for line in f:
            g = line.split()[0]
            if g not in genes:
                genes.append(g)
    gene_list.append(genes)

header = ["Target", "Accession",] + labels
print "\t".join(header)
for g in gene_list[0]:
    row = g.split('.') + ['Y']
    for genes in gene_list[1:]:
        if g in genes:
            row.append('Y')
        else:
            row.append('N')
    print "\t".join(row)