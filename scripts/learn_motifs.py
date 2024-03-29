#!/usr/bin/env python

# Imports
import numpy as np
import sys

try:
    import cPickle as pickle
except:
    import pickle

from strum import strum

# Set up indexes for defining the weight matrices
nucs = "ACGT"
nuc_index = dict(zip(nucs, range(4)))
dimers = [a+b for a in nucs for b in nucs]
di_index = dict(zip(dimers, range(16)))
N_seq = 500

def learn(basename, n_process, random_seed):
    print >> sys.stderr, "Read in MEME seqs"
    meme_output = open("output/" + basename + "_meme/meme_out/meme.txt")

    # Load in the sequences determined by MEME as sufficient to train the PWM
    start = False
    count = 0
    sequences = []
    for line in meme_output:
        if "1 sites sorted by position p-value" in line:
            start = True
        if start:
            count += 1
            if count >= 5:
                if "--------------" in line:
                    break
                sequence = line.split()[5]
                sequences.append(sequence.upper())

    # Train PWM, DWM, and ML-StruM on sequences from MEME output
    print >> sys.stderr, "Train ML motifs"
    pwm = learn_pwm(sequences)
    dwm = learn_dwm(sequences)
    ml_strum = learn_strum(sequences, n_process)
    ml_strum.filter()

    k = pwm.shape[1]

    # Load in original sequences from ChIP experiment
    sequences = fasta_reader(open("data/{}.fa".format(basename)))

    # Train EM-StruM on first N_seq peaks from ChIP experiment
    print >> sys.stderr, "Train EM StruM"
    prelim_em_strum = strum.StruM(mode='proteingroove', n_process=n_process)
    prelim_em_strum.train_EM(sequences[:N_seq], fasta=False, lim=0.1, 
        k=(k-1), max_iter=250, n_init=30, random_seed=random_seed,)
    prelim_kmers = []
    for seq in sequences[:N_seq]:
        rseq = strum.rev_comp(seq)
        s1 = prelim_em_strum.score_seq(seq)
        s2 = prelim_em_strum.score_seq(rseq)
        if np.max(s1) > np.max(s2):
            i = np.argmax(s1)
        else:
            seq = rseq
            i = np.argmax(s2)
        prelim_kmers.append(seq[i:i+k])

    em_strum = strum.StruM(mode='full', n_process=n_process)
    em_strum.train(prelim_kmers, lim=0.1)
    em_strum.filter()

    return (pwm, dwm, ml_strum, em_strum), sequences

def learn_pwm(sequences):
    """Generate a PWM from a set of training sequences."""
    """PWM = Position weight matrix"""
    PWM = np.zeros([4, len(sequences[0])])
    for seq in sequences:
        for i,n in enumerate(seq):
            PWM[nuc_index[n], i] += 1
    PWM += 1
    PWM /= (len(sequences) + 4.)
    return PWM

def learn_dwm(sequences):
    """Generate a DWM from a set of training sequences."""
    """DWM = Dinucleotide weight matrix"""
    DWM = np.zeros([16, len(sequences[0])])
    for seq in sequences:
        for i,n in enumerate(seq):
            if i == 0:
                DWM[nuc_index[n], i] += 1
            else:
                di = seq[i-1:i+1]
                DWM[di_index[di], i] += 1
    DWM += 1
    DWM[:,0] /= (len(sequences) + 4.)
    DWM[:, 1:] /= (len(sequences) + 16.)
    return DWM

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

def learn_strum(sequences, n_process):
    """Generate a StruM from a set of training sequences."""
    """StruM = Structural Motif"""
    motif = strum.StruM(mode='full', n_process=n_process)
    motif.train(sequences, lim=0.1)
    return motif

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print "Usage:\n\t{} <basename> <n_process> <random_seed> <output.p>"
        quit()

    basename = sys.argv[1]
    n_process = int(sys.argv[2])
    random_seed = int(sys.argv[3])
    models, sequences = learn(basename, n_process, random_seed)
    pickle.dump(models, open(sys.argv[4], 'wb'))
