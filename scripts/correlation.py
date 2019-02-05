#!/usr/bin/env python

import numpy as np
import random
import sys

#from sklearn.linear_model import LogisticRegression as logit

try:
	import cPickle as pickle
except:
	import pickle

# import StruM
N = 1000

# Set up indexes for defining the weight matrices
nucs = "ACGT"
nuc_index = dict(zip(nucs, range(4)))
dimers = [a+b for a in nucs for b in nucs]
di_index = dict(zip(dimers, range(16)))

def main(tf):
	pwm, dwm, ml_strum, em_strum, logit, logit2 = pickle.load(open("output/{}.p".format(tf), "rb"))

	k = pwm.shape[1]
	sequence = "".join([random.choice(nucs) for i in range(N+k)])
	sequences = [sequence[i:i+k] for i in range(N)]

	scores = []
	for kmer in sequences:
		s1 = score_PWM(pwm, kmer)
		s2 = score_DWM(dwm, kmer)
		scores.append((s1, s2))
	s3 = ml_strum.score_seq_filt(sequence)
	s4 = em_strum.score_seq_filt(sequence)
	

	scores = np.vstack(scores).T
	scores = np.vstack([scores, np.reshape(s3, [1,-1]), np.reshape(s4, [1,-1])])
	coef = np.corrcoef(scores)

	idx = np.triu_indices(4, 1)
	print "{}\t".format(tf) + "\t".join(["{:0.3f}".format(x) for x in coef[idx]])

def score_PWM(PWM, kmer):
	"""Score a kmer with a given PWM, and return log2 of the score."""
	p = np.sum([ PWM[nuc_index[n],j] for j,n in enumerate(kmer)])
	return p
	#return np.log(p/np.product([0.25]*PWM.shape[1]))
	return np.log2(p)

def score_DWM(DWM, kmer):
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

def score_StruM(strum, kmer):
	"""Score a kmer with a given StruM, and return log10 of the score."""
	return strum.eval_filt(strum.translate(kmer))

if __name__ == '__main__':
	tf = sys.argv[1]
	main(tf)