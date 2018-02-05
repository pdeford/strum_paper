#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import random
import sys
from sklearn.metrics import roc_curve, roc_auc_score

from strum import strum
import bx.bbi.bigwig_file


DNase_bigwig_path = sys.argv[1]
chrom_path = sys.argv[2]
factor = sys.argv[3]

training_bed = open("output/unique_{}.bed".format(factor))
pos_bed = open("output/med_{}.bed".format(factor))
neg_bed = open("output/not_K562_{}.bed".format(factor))

DNASE_TRAIN  = open("output/train_dnase_{}.txt".format(factor),  "wb")
DNASE_TEST   = open("output/test_dnase_{}.txt".format(factor),   "wb")
LABEL_TRAIN  = open("output/train_label_{}.txt".format(factor),  "wb")
LABEL_TEST   = open("output/test_label_{}.txt".format(factor),   "wb")
MOTIF_TRAIN  = open("output/train_motif_{}.txt".format(factor),  "wb")
MOTIF_TEST   = open("output/test_motif_{}.txt".format(factor),   "wb")
MOTIF_TRAIN2 = open("output/train_motif2_{}.txt".format(factor), "wb")
MOTIF_TEST2  = open("output/test_motif2_{}.txt".format(factor),  "wb")

PWM_OUTPUT_TRAIN = open("output/train_pwm_stuff_{}.txt".format(factor), "wb")
PWM_OUTPUT_TEST  = open("output/test_pwm_stuff_{}.txt".format(factor),  "wb")

PWM_DNASE_TRAIN  = open("output/train_pwm_dnase_{}.txt".format(factor), "wb")
PWM_DNASE_TEST   = open("output/test_pwm_dnase_{}.txt".format(factor),  "wb")

AUC1_OUT = open("output/auc1_{}.txt".format(factor), "wb")
AUC2_OUT = open("output/auc2_{}.txt".format(factor), "wb")
AUC3_OUT = open("output/auc3_{}.txt".format(factor), "wb")

def lookup_sequence(chrom,start=None,end=None):
    if start and end:
        with open(chrom_path+chrom+".fa") as f_c:
            f_c.readline()
            offset = f_c.tell()
            f_c.seek(start + offset + (start//50))
            lines = f_c.read(end-start + (end-start)//50 + 1)
        sequence = "".join(lines.split("\n"))
    else:
        with open(chrom_path+chrom+".fa") as f_c:
            sequence = "".join(f_c.read().split("\n")[1:])
    return sequence.upper()

def lookup_DNase(data, chrom, start, end, extend = False):
	bwh = bx.bbi.bigwig_file.BigWigFile(open(data))
	if extend:
		extend = abs(start-end)-1
	else:
		extend = 0
	trace = bwh.get_as_array(chrom, min(start,end)-extend, max(start, end)-1+extend)
	if trace is None:
		trace = np.zeros(abs(start-end)-1+2*extend)
	trace[np.isnan(trace)] = 0.0
	trace -= np.min(trace)

	if start > end:
		trace = trace[::-1]

	if not all(trace==0.0):
		trace /= np.max(trace)

	if extend != 0:
		return np.reshape(trace, [3,-1])
	else:
		return trace

def bed2seq(bedfile, n_sequences=200, pval_col=6):
	seqs = []
	scores = []
	positions = []
	for line in bedfile:
		fields = line.split()
		chrom, start, end = fields[0], int(fields[1]), int(fields[2])
		s = lookup_sequence(chrom, start, end)
		seqs.append(s)
		scores.append(float(fields[pval_col]))
		positions.append((chrom, start, end))
	return ([ x[0] for x in sorted(zip(seqs, scores), key=lambda x:x[1], reverse=True) ][:min(len(scores), n_sequences)],
	        [ x[0] for x in sorted(zip(positions, scores), key=lambda x:x[1], reverse=True) ][:min(len(scores), n_sequences)])

##################################################
# LOAD PWM FROM MEME OUTPUT
f_meme = open("output/{}_meme/meme_out/meme.txt".format(factor))
PWM = []
start = False
count = -3
for row in f_meme:
	if "Motif 1 position-specific probability matrix" in row:
		start = True
	if start:
		count += 1
		if count > 0:
			if "------" in row:
				break
			pos = [ float(x) for x in row.split() ]
			PWM.append(pos)
PWM = np.asarray(PWM).T
##################################################

nuc_index = dict(zip("ACGT", range(4)))

def score_seq_p(PWM, seq):
	k = PWM.shape[1]
	f_scores = []
	r_scores = []
	for i in range(len(seq) - k + 1):
		f, r = 1, 1
		kmer = seq[i:i+k]
		rkmer = motif.rev_comp(kmer)
		for j,n in enumerate(kmer):
			f *= PWM[nuc_index[n], j]
			r *= PWM[nuc_index[rkmer[j]], j]
		f_scores.append(f)
		r_scores.append(r)
	return f_scores, r_scores


print "Extract sequences"
little_n = 200
training_sequences, training_positions = bed2seq(bedfile=training_bed, n_sequences=little_n, pval_col=7)
pos_sequences, pos_positions = bed2seq(bedfile=pos_bed, n_sequences=2000, pval_col=7)
neg_sequences, neg_positions = bed2seq(bedfile=neg_bed, n_sequences=2000, pval_col=7)

test_sequences = pos_sequences + neg_sequences
test_positions = pos_positions + neg_positions
test_sequences2 = training_sequences + neg_sequences[-little_n:]
test_positions2 = training_positions + neg_positions[-little_n:]

Y  = [1 for x in pos_sequences] + [0 for x in neg_sequences]
Y2 = [1 for x in training_sequences] + [0 for x in neg_sequences[:little_n]]


print "Train StruM"
motif = strum.StruM(load_diprodb=True, mode="protein", n_process=-1)
motif.train_EM(training_sequences, fasta=False, k=10, max_iter=200, random_seed=808, lim=0.001, n_init=5)

motif.print_PWM(labels=True)


print "Find best scoring position in each training peak"
positions = []
weights = []
strand = []
pwm_scores = []
pwm_positions = []
pwm_strand = []

for seq in training_sequences:
	rseq = motif.rev_comp(seq)
	f_scores = motif.score_seq(seq)
	r_scores = motif.score_seq(rseq)
	f_i = np.argmax(f_scores)
	r_i = np.argmax(r_scores)
	if f_scores[f_i] > r_scores[r_i]:
		positions.append(f_i)
		weights.append(f_scores[f_i])
		strand.append(1)
	else:
		positions.append(len(seq)-r_i-motif.k)
		weights.append(r_scores[r_i])
		strand.append(-1)
	f_scores, r_scores = score_seq_p(PWM, seq)
	f_i = np.argmax(f_scores)
	r_i = np.argmax(r_scores)
	if f_scores[f_i] > r_scores[r_i]:
		pwm_positions.append(f_i)
		pwm_scores.append(f_scores[f_i])
		pwm_strand.append(1)
	else:
		pwm_positions.append(len(seq)-r_i-motif.k)
		pwm_scores.append(r_scores[r_i])
		pwm_strand.append(-1)
MOTIF_TRAIN.write(' '.join([str(x) for x in weights]) + "\n")
LABEL_TRAIN.write(' '.join([str(1) for x in weights]) + "\n")
PWM_OUTPUT_TRAIN.write(' '.join([str(x) for x in pwm_scores]) + "\n")


print "Get standard scores for each test region"
X1 = []
positions_X1 = []
strand_X1 = []
pwm_scores_X1 = []
pwm_positions_X1 = []
pwm_strand_X1 = []
for seq in test_sequences:
	rseq = motif.rev_comp(seq)
	f_scores = motif.score_seq(seq)
	r_scores = motif.score_seq(rseq)
	X1.append(np.max(np.hstack([f_scores,r_scores])))
	f_i = np.argmax(f_scores)
	r_i = np.argmax(r_scores)
	if f_scores[f_i] > r_scores[r_i]:
		positions_X1.append(f_i)
		#X1.append(f_scores[f_i])
		strand_X1.append(1)
	else:
		positions_X1.append(len(seq)-r_i-motif.k)
		#X1.append(r_scores[r_i])
		strand_X1.append(-1)
	f_scores, r_scores = score_seq_p(PWM, seq)
	f_i = np.argmax(f_scores)
	r_i = np.argmax(r_scores)
	if f_scores[f_i] > r_scores[r_i]:
		pwm_positions_X1.append(f_i)
		pwm_scores_X1.append(f_scores[f_i])
		pwm_strand_X1.append(1)
	else:
		pwm_positions_X1.append(len(seq)-r_i-motif.k)
		pwm_scores_X1.append(r_scores[r_i])
		pwm_strand_X1.append(-1)
MOTIF_TEST.write(' '.join([str(x) for x in X1]) + "\n")
LABEL_TEST.write(' '.join([str(x) for x in Y]) + "\n")
PWM_OUTPUT_TEST.write(' '.join([str(x) for x in pwm_scores_X1]) + "\n")


pwm_k = PWM.shape[1]
print "Extracting DNase signal"
DNase_signals = []
for i, (chrom, start, stop) in enumerate(training_positions):
	addition = positions[i]
	new_start = start + addition
	trace = lookup_DNase(DNase_bigwig_path, chrom, new_start, new_start + motif.k, True ).ravel()
	if strand[i] == -1:
		trace = trace[::-1]
	DNase_signals.append(trace)
	if strand[i] == 1:
		DNASE_TRAIN.write(" ".join([str(x) for x in lookup_DNase(DNase_bigwig_path, chrom, new_start-50, new_start+motif.k+50)]) + "\n")
	else:
		DNASE_TRAIN.write(" ".join([str(x) for x in lookup_DNase(DNase_bigwig_path, chrom, new_start+motif.k+50, new_start-50)]) + "\n")
	
	addition = pwm_positions[i]
	new_start = start + addition
	if pwm_strand[i] == 1:
		PWM_DNASE_TRAIN.write(" ".join([str(x) for x in lookup_DNase(DNase_bigwig_path, chrom, new_start-50, new_start+pwm_k+50)]) + "\n")
	else:
		PWM_DNASE_TRAIN.write(" ".join([str(x) for x in lookup_DNase(DNase_bigwig_path, chrom, new_start+pwm_k+50, new_start-50)]) + "\n")



print "Learn DNase component of StruM"
print motif.k
DNase_signals = np.asarray(DNase_signals)[:, :3*motif.k]
strum_addition = [np.reshape(np.average(DNase_signals, axis=0), [3,-1]), np.reshape(np.std(DNase_signals, axis=0), [3,-1])]



print "Add DNase to StruM"
avgs = np.hstack([np.reshape(motif.strum[0], [1, -1]),
				  np.reshape(strum_addition[0], [1, -1])
				  ])
stds = np.hstack([np.reshape(motif.strum[1], [1, -1]),
				  np.reshape(strum_addition[1], [1, -1])
				  ])
stds[stds < 0.001] = 0.001

motif.strum = [avgs, stds]
motif.update(data=DNase_bigwig_path, func=lookup_DNase, features=['DNaseUpstream', 'DNaseCenter', 'DNaseDownstream'])


print "Get DNase scores for each test region"
X2 = []
for i, seq in enumerate(test_sequences):
	chrom, start, stop = test_positions[i]
	seq = seq[:stop-start]
	rseq = motif.rev_comp(seq)
	f_scores = motif.score_seq(seq,  *(chrom, start, stop, True))
	r_scores = motif.score_seq(rseq, *(chrom, stop, start, True))
	X2.append(np.max(np.hstack([ f_scores, r_scores ])))

	addition = positions_X1[i]
	new_start = start + addition
	if strand_X1[i] == 1:
		DNASE_TEST.write(" ".join([str(x) for x in lookup_DNase(DNase_bigwig_path, chrom, new_start-50, new_start+motif.k+50)]) + "\n")
	else:
		DNASE_TEST.write(" ".join([str(x) for x in lookup_DNase(DNase_bigwig_path, chrom, new_start+motif.k+50, new_start-50)]) + "\n")

	addition = pwm_positions_X1[i]
	new_start = start + addition
	if pwm_strand_X1[i] == 1:
		PWM_DNASE_TEST.write(" ".join([str(x) for x in lookup_DNase(DNase_bigwig_path, chrom, new_start-50, new_start+pwm_k+50)]) + "\n")
	else:
		PWM_DNASE_TEST.write(" ".join([str(x) for x in lookup_DNase(DNase_bigwig_path, chrom, new_start+pwm_k+50, new_start-50)]) + "\n")

X3 = []
for i, seq in enumerate(test_sequences2):
	chrom, start, stop = test_positions2[i]
	seq = seq[:stop-start]
	rseq = motif.rev_comp(seq)
	f_scores = motif.score_seq(seq,  *(chrom, start, stop, True))
	r_scores = motif.score_seq(rseq, *(chrom, stop, start, True))
	X3.append(np.max(np.hstack([ f_scores, r_scores ])))

X1, X2, X3 = np.asarray(X1), np.asarray(X2), np.asarray(X3)
X4 = np.asarray(pwm_scores_X1)


print "ROC Curve creation"
X1 = X1 - np.min(X1)
X1 = X1 / np.max(X1)

X2 = X2 - np.min(X2)
X2 = X2 / np.max(X2)

X3 = X3 - np.min(X3)
X3 = X3 / np.max(X3)

X4 = X4 - np.min(X4)
X4 = X4 / np.max(X4)

# Make ROC Curve
auc1 = roc_auc_score(Y, X1)
fpr1, tpr1, thresh = roc_curve(Y, X1)
print "auROC curve no DNase: {} ".format(auc1)

auc2 = roc_auc_score(Y, X2)
fpr2, tpr2, thresh = roc_curve(Y, X2)
print "auROC curve w/ DNase: {} ".format(auc2)

auc3 = roc_auc_score(Y2, X3)
fpr3, tpr3, thresh = roc_curve(Y2, X3)
print "auROC curve w/ DNase and training sequences: {} ".format(auc3)

auc4 = roc_auc_score(Y, X4)
fpr4, tpr4, thresh = roc_curve(Y, X4)
print "auROC curve w/ PWM: {}".format(auc4)

# plt.figure( figsize=[10, 10] )
# plt.plot( [0, 1], [0, 1], '--', color='gray' )
# plt.plot( fpr1, tpr1, label="StruM (AUC = {:0.2f})".format(auc1) )
# plt.plot( fpr2, tpr2, label="DNase-StruM (AUC = {:0.2f})".format(auc2) )
# plt.plot( fpr3, tpr3, label="DNase-StruM control (AUC = {:0.2f})".format(auc3) )
# plt.plot( fpr4, tpr4, label="PWM (AUC = {:0.2f})".format(auc4) )
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.legend(loc='lower right')
# plt.savefig('output/DNase_performance.png')
# plt.close()

AUC1_OUT.write(" ".join([str(x) for x in fpr1]) + "\n")
AUC2_OUT.write(" ".join([str(x) for x in fpr2]) + "\n")
AUC3_OUT.write(" ".join([str(x) for x in fpr4]) + "\n")
AUC1_OUT.write(" ".join([str(x) for x in tpr1]) + "\n")
AUC2_OUT.write(" ".join([str(x) for x in tpr2]) + "\n")
AUC3_OUT.write(" ".join([str(x) for x in tpr4]) + "\n")
