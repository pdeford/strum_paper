#!/usr/bin/env python

# Imports
import matplotlib as mpl
mpl.use('Agg') # Allows it to run on the server

import matplotlib.pyplot as plt
import numpy as np
import random
import sys

from multiprocessing import Pool
from scipy import interp
from sklearn.linear_model import LogisticRegression as logit
from sklearn.metrics import roc_curve, auc
from sklearn.cross_validation import StratifiedKFold

try:
	import cPickle as pickle
except:
	import pickle

from strum import strum

def score_all(basename, n_process, random_seed, models):
	global k, pwm, dwm, ml_strum, em_strum

	random.seed(random_seed)

	pwm, dwm, ml_strum, em_strum = models
	k = pwm.shape[1]

	print >> sys.stderr, "Load positive and negative sequences"
	# Get positive examples from file
	peaks = sequences[N_seq:]

	# Create background sequences for classification from ChIP sequences
	decoys = [list(p) for p in peaks]
	for p in decoys:
		random.shuffle(p)
	decoys = ["".join(p) for p in decoys]

	print >> sys.stderr, "Score with 4 motifs"
	# Score each sequence with each motif
	y = []
	data = []
	for i, seq_set in enumerate([decoys, peaks]):
		y += [i] * len(seq_set)
		pool = Pool(n_process)
		data_ = pool.map(eval_seq, seq_set)
		pool.close()
		pool.join()
		data += data_

	data = np.vstack(data)
	y = np.array(y, dtype=int)

	print >> sys.stderr, "Get ROC curves"
	# Evaluate the performance of each motif
	results = []
	for i, motif in enumerate([pwm, dwm, ml_strum, em_strum]):
		fpr, tpr, _ = roc_curve(y, data[:,i])
		results.append((fpr, tpr, auc(fpr,tpr)))

	############## REPEAT FOR OTHER TYPES OF DECOYS #############
	print >> sys.stderr, "Load other decoys"
	decoys2 = fasta_reader(open("data/{}.flank.fa".format(basename)))[N_seq:]

	print >> sys.stderr, "Score decoys"
	pool = Pool(n_process)
	data2 = pool.map(eval_seq, decoys2)
	pool.close()
	pool.join()
	data2 = np.vstack(data2)

	print >> sys.stderr, "Decoys ROCs"
	results2 = []
	working_x = np.vstack([data[y == 1], data2])
	working_y = np.zeros([len(working_x),], dtype=int)
	working_y[:sum(y)] = 1
	for i, motif in enumerate([pwm, dwm, ml_strum, em_strum]):
		fpr, tpr, _ = roc_curve(working_y, working_x[:,i])
		results2.append((fpr, tpr, auc(fpr,tpr)))

	##############################################################

	print >> sys.stderr, "Logit"
	# Train a combined model using logistic regression
	## Scale the input so each set of scores ranges in [0,1]
	scaler = (np.min(data, axis=0), np.max(data, axis=0))
	data -= scaler[0]
	data /= (scaler[1] - scaler[0])

	## Train the model
	clf = logit()
	cv = StratifiedKFold(y,n_folds=10)
	tprs = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100)
	for train, test in cv:
		probas_ = clf.fit(data[train], y[train]).predict_proba(data[test])
		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
		tprs.append(interp(mean_fpr, fpr, tpr))
		tprs[-1][0] = 0.0
		roc_auc = auc(fpr, tpr)
		aucs.append(roc_auc)
	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	results.append((mean_fpr, mean_tpr, np.mean(aucs))) # Combined model

	############## REPEAT FOR OTHER TYPES OF DECOYS #############
	# Train a combined model using logistic regression
	## Scale the input so each set of scores ranges in [0,1]
	print >> sys.stderr, "...and Again"
	scaler2 = (np.min(working_x, axis=0), np.max(working_x, axis=0))
	working_x -= scaler2[0]
	working_x /= (scaler2[1] - scaler2[0])

	## Train the model
	clf2 = logit()
	cv2 = StratifiedKFold(working_y, n_folds=10)
	tprs = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100)
	for train, test in cv2:
		probas_ = clf2.fit(working_x[train], working_y[train]).predict_proba(working_x[test])
		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(working_y[test], probas_[:, 1])
		tprs.append(interp(mean_fpr, fpr, tpr))
		tprs[-1][0] = 0.0
		roc_auc = auc(fpr, tpr)
		aucs.append(roc_auc)
	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	results2.append((mean_fpr, mean_tpr, np.mean(aucs))) # Combined model
	##############################################################

	print >> sys.stderr, "Print Results"
	# Output results and plot the results
	print basename + "\t" + "\t".join(["%0.6f" % x[-1] for x in results])
	############## REPEAT FOR OTHER TYPES OF DECOYS #############
	print basename + "\t" + "\t".join(["%0.6f" % x[-1] for x in results2])
	#############################################################

	labels = ["PWM",       "DWM",        "ML-StruM", "EM-StruM",  "Combined\nLOGIT"]
	colors = ["firebrick", "darkorange", "seagreen", "steelblue", "mediumpurple"]

	print >> sys.stderr, "Generate Plots"
	plt.figure()
	for i, name in enumerate(labels):
		row = results[i]
		plt.plot(row[0], row[1], c=colors[i], label=name)
	plt.xlim([-0.05, 1.05])
	plt.ylim([-0.05, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	#plt.title('Receiver operating characteristic example')
	plt.legend(loc="lower right")
	plt.savefig("output/{}_ROC.pdf".format(basename))
	plt.close()

	plt.figure()
	for i, name in enumerate(labels):
		row = results[i]
		plt.plot(row[0], row[1], c=colors[i], label=name)
	plt.xlim([-0.05, 1.05])
	plt.ylim([-0.05, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	#plt.title('Receiver operating characteristic example')
	plt.legend(loc="lower right")
	plt.savefig("output/{}_ROC2.pdf".format(basename))
	plt.close()

	return (clf, scaler), (clf2, scaler2), results, results2


def score_pwm(PWM, kmer):
	"""Score a kmer with a given PWM, and return log2 of the score."""
	p = np.sum([ PWM[nuc_index[n],j] for j,n in enumerate(kmer)])
	return p
	#return np.log(p/np.product([0.25]*PWM.shape[1]))
	return np.log2(p)

def score_dwm(DWM, kmer):
	"""Score a kmer with a given DWM, and return log2 of the score."""
	p = DWM[nuc_index[kmer[0]],0]
	for j,n in enumerate(kmer):
		if j == 0:
			continue
		di = kmer[j-1:j+1]
		p += DWM[di_index[di],j]
	return p
	#return np.log(p/np.product([0.25]*PWM.shape[1]))
	return np.log2(p)

def score_strum(strum, kmer):
	"""Score a kmer with a given StruM, and return log10 of the score."""
	return strum.eval(strum.translate(kmer))

def rev_comp(seq):
	"""Return the reverse complement of a nucleotide sequence."""
	nucs = "ACNGT"
	index = dict(zip(nucs, nucs[::-1]))
	return "".join([index[n] for n in seq][::-1])

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

def eval_seq(seq):
		pwm_scores = []
		dwm_scores = []
		ml_strum_scores = []
		em_strum_scores = []
		for j in range(len(seq) - k + 1):
			kmer = seq[j:j+k]
			if 'N' in kmer: continue
			rkmer = rev_comp(kmer)
			pwm_scores.append(score_pwm(pwm, kmer))
			pwm_scores.append(score_pwm(pwm, rkmer))
			dwm_scores.append(score_dwm(dwm, kmer))
			dwm_scores.append(score_dwm(dwm, rkmer))
			ml_strum_scores.append(score_strum(ml_strum, kmer))
			ml_strum_scores.append(score_strum(ml_strum, rkmer))
			em_strum_scores.append(score_strum(em_strum, kmer))
			em_strum_scores.append(score_strum(em_strum, rkmer))
		return (np.max(pwm_scores), np.max(dwm_scores), np.max(ml_strum_scores), np.max(em_strum_scores))

if __name__ == '__main__':
	basename = sys.argv[1]
	n_process = int(sys.argv[2])
	random_seed = int(sys.argv[3])
	models = pickle.load(open(sys.argv[4], 'rb'))
	
	score_all(basename, n_process, random_seed, models)