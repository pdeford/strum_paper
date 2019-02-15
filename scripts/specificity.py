import sys
import numpy as np
from sklearn.metrics import roc_curve, auc, precision_recall_curve

try:
    import cPickle as pickle
except:
    import pickle

basename = sys.argv[1]
seed = int(sys.argv[2])
n_process= int(sys.argv[3])

tf, accession = basename.split('.')
np.random.seed(seed)

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

def rev_comp(seq):
    """Return the reverse complement of a nucleotide sequence."""
    nucs = "ACNGT"
    index = dict(zip(nucs, nucs[::-1]))
    return "".join([index[n] for n in seq][::-1])

f = open("output/families.txt")

families = []
for line in f:
    fields = line.split('\t')
    family = fields[1]
    if family not in families:
        families.append(family)
    if tf == fields[0]:
        tf_fam = fields[1]


try:
    with open('data/{}.fa'.format(basename)) as f:
        positives = fasta_reader(f)[500:1000]
except:
    quit()

sequences = []
for family in sorted(families):
    if family != tf_fam:
        try:
            with open("data/{}.fa".format(family)) as f:
                sequences += fasta_reader(f)
        except:
            continue

negatives = list(np.random.choice(sequences, 500, False))

with open("data/{}.fa".format(tf_fam)) as f:
    sequences = fasta_reader(f)

negatives2 = list(np.random.choice(sequences, 500, False))


pwm, dwm, ml_strum, em_strum, (logit, scaler), (logit2, scaler2) = pickle.load(open("output/{}.p".format(basename), "rb"))

y = []
y2 = []
x = []
x2 = []
for i, seq_set in enumerate([negatives, positives]):
    for seq in seq_set:
        y.append(i)
        rseq = rev_comp(seq)
        s1 = np.max(em_strum.score_seq_filt(seq))
        s2 = np.max(em_strum.score_seq_filt(rseq))
        x.append(max(s1,s2))
        if i == 1:
            y2.append(i)
            x2.append(max(s1,s2))

for seq in negatives2:
    rseq = rev_comp(seq)
    s1 = np.max(em_strum.score_seq_filt(seq))
    s2 = np.max(em_strum.score_seq_filt(rseq))
    x2.append(max(s1,s2))    
    y2.append(0)

fpr, tpr, _ = roc_curve(y, x)
precision, recall, _ = precision_recall_curve(y, x)
results = (fpr, tpr, auc(fpr,tpr), recall, precision, auc(recall, precision))

fpr2, tpr2, _ = roc_curve(y2, x2)
precision2, recall2, _ = precision_recall_curve(y2, x2)
results2 = (fpr2, tpr2, auc(fpr2,tpr2), recall2, precision2, auc(recall2, precision2))

print basename, results[2], results[5], results2[2], results2[5]