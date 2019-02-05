#!/usr/bin/env python

# Imports
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import sys
import subprocess
# plt.style.use("dark_background")

from scipy.stats import pointbiserialr

try:
        import cPickle as pickle
except:
        import pickle

# Import StruM
from strum import strum

# Prepare data
nucs = "ACGT"
nuc_index = dict(zip(nucs, range(4)))
dimers = [a+b for a in nucs for b in nucs]
di_index = dict(zip(dimers, range(16)))

# Define Main function
def main(basename, n_process=1):
	global k, pwm, em_strum

	# Load precomputed TFBS models
	models = pickle.load(open("output/{}.p".format(basename)))
	pwm = models[0]
	em_strum = models[3]

	k = pwm.shape[1]

	f = open("output/{}_meme/fimo_out_1/fimo.txt".format(basename))
	f.readline()
	matches = {}
	for line in f:
		fields = line.split()
		chrom = fields[2]
		start = int(fields[3])
		matches.setdefault(chrom, []).append(start)

	if sum([len(matches[chrom]) for chrom in matches]) < 100:
		print >> sys.stderr, "Not enough significance (n < 100). Exiting..."
		return

	for chrom in matches.keys():
		matches[chrom].sort()

	regions = {}
	for chrom in matches:
		matchlist = matches[chrom]
		tmp = [matchlist[0]]
		for i in range(1, len(matchlist)):
			if matchlist[i] - matchlist[i-1] < 10:
				tmp.append(matchlist[i])
			else:
				regions.setdefault(chrom, []).append(int(np.average(tmp)))
				tmp = [matchlist[i]]
		else:
			regions.setdefault(chrom, []).append(int(np.average(tmp)))

	with open("output/{}_top_matches.bed".format(basename), "wb") as f:
		for chrom in sorted(regions.keys()):
			for pos in regions[chrom]:
				start = pos - 100
				if start < 0:
					continue
				f.write("{}\t{}\t{}\n".format(chrom, start, pos+100+k))

	subprocess.call(
		"bedtools getfasta -fi data/hg19.fa -bed {} -fo 'data/{}_matches.fa'".format(
			"output/{}_top_matches.bed".format(basename), basename), 
		shell=True)

	# Load in original sequences from ChIP experiment
	headers, sequences = fasta_reader(open("data/{}_matches.fa".format(basename)))
	n = sum([len(matches[chrom]) for chrom in matches])

	####SCORE THE SEQUENCES WITH STRUM
	####FILTER FOR `n` TOP SITES
	####LOOK AT POSITIONS RELATIVE TO FIMO POSITIONS
	scores = []
	for seq in sequences:
		scores.append(eval_seq(seq))
	scores = np.vstack(scores)

	# BINARIZE THE STRUM SCORES
	for thresh in np.sort(scores.ravel())[::-1]:
		count = np.sum(scores > thresh)
		if count > n:
			break
	
	bin_scores = np.zeros(scores.shape)
	bin_scores[scores > thresh] = 1

	fimo_mat = np.zeros(scores.shape)
	adjust = 0
	i = 0
	for chrom in sorted(regions.keys()):
		for r in regions[chrom]:
			left = r - 100
			if left < 0:
				continue
			right = r + 100
			for pos in matches[chrom]:
				if pos < left:
					continue
				elif pos >= right: 
					break
				else:
					j = pos - left
					fimo_mat[i, j] = 1
			i += 1

	# SORT THE MATRICES FOR THE HEATMAP
	corrs = []
	for i in range(fimo_mat.shape[0]):
		r,p = pointbiserialr(fimo_mat[i], scores[i])
		corrs.append(r)
	sort_idx = np.argsort(corrs)
	fimo_mat = fimo_mat[sort_idx]
	scores = scores[sort_idx]
	bin_scores = bin_scores[sort_idx]

	# for i,r in enumerate(regions):
	# 	left = r - 100
	# 	if left < 0: 
	# 		adjust += 1
	# 		continue
	# 	right = r + 100
	# 	for pos in matches:
	# 		if pos < left:
	# 			continue
	# 		elif pos >= right:
	# 			break
	# 		else:
	# 			j = pos - left
	# 			fimo_mat[i-adjust,j] = 1

	plt.figure(figsize=[12,6])
	plt.subplot(121)
	plt.pcolor(fimo_mat, cmap='YlGnBu')
	plt.axis('tight')
	# plt.xticks([0, 49, 99, 149, 199], [-100, -50, 0, 50, 100])
	# plt.xlabel("Position (bp)")
	plt.yticks([])
	plt.xticks([])
	plt.subplot(122)
	#plt.pcolor(bin_scores, cmap='YlGnBu')
	plt.pcolor(scores, cmap='YlGnBu')
	plt.axis('tight')
	# plt.xticks([0, 49, 99, 149, 199], [-100, -50, 0, 50, 100])
	# plt.xlabel("Position (bp)")
	plt.yticks([])
	plt.xticks([])
	plt.subplots_adjust(left=0,bottom=0,top=1,right=1,wspace=0.01)
	plt.savefig("output/{}_fimo_v_strum_matches.png".format(basename))
	plt.close()

	diffs = []
	for i, row in enumerate(bin_scores):
		row2 = fimo_mat[i]
		for m1 in np.where(row == 1)[0]:
			delta = np.absolute(m1 - np.where(row2 == 1)[0])
			diffs.append(min(delta))

	plt.figure()
	plt.hist(diffs, bins=30)
	ax = plt.gca()
	ax.tick_params(labelleft='off') 
	plt.xlabel("Distance to nearest match (bp)")
	plt.savefig("output/{}_fimo_v_strum_dist.pdf".format(basename))
	plt.close()

	print basename, n, k, np.average(diffs), np.std(diffs)


# Define Functions
def rev_comp(seq):
	"""Return the reverse complement of a nucleotide sequence."""
	nucs = "ACNGT"
	index = dict(zip(nucs, nucs[::-1]))
	return "".join([index[n] for n in seq][::-1])

def fasta_reader(file_obj):
	"""Read a FASTA file and return the sequences in a list."""
	headers = []
	sequences = []
	currSeq = ""
	firstLine = True
	for line in file_obj:
		if line.startswith(">") or line.strip() == "":
			header = line.strip()[1:]
			if firstLine:
				firstLine = False
			else:
				sequences.append(currSeq.upper())
				currSeq = ""
			headers.append(header)
		else:
			s = line.strip()
			currSeq += s
		if line == "":
			break
	if currSeq != "":
		sequences.append(currSeq.upper())
	return headers[:-1], sequences

def eval_seq(seq):
	em_strum_scores_f = em_strum.score_seq_filt(seq)
	em_strum_scores_r = em_strum.score_seq_filt(rev_comp(seq))[::-1]
	em_strum_scores = np.max(np.vstack([em_strum_scores_f, em_strum_scores_r]), axis=0)
	return em_strum_scores

if __name__ == '__main__':
	basename = sys.argv[1]
	n_process = int(sys.argv[2])
	main(basename, n_process)
