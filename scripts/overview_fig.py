#!/usr/bin/env python

import matplotlib as mpl
mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

pwm_c = """   1   2   3   4   5   6   7   8   9
A 0.2 0.7 0.0 1.0 0.0 1.0 0.9 0.2 0.4
C 0.4 0.0 0.0 0.0 0.0 0.0 0.0 0.2 0.2
G 0.3 0.0 1.0 0.0 0.0 0.0 0.1 0.6 0.3
T 0.1 0.3 0.0 0.0 1.0 0.0 0.0 0.0 0.1"""

def cm2inch(value):
    return value/2.54

onecol = cm2inch(84/10.)
twocol = cm2inch(178/10.)
maxheight = cm2inch(230/10.)

def plotter():
	x_vals = np.linspace(0, 28, 500)
	helix_a = np.sin(x_vals)
	helix_b = np.sin(x_vals + 2)
	sequences = [
		"CTGACAGATAAGACTATGACGGTACC",
		"TGCTGATCGGTACGAGATAAGAGGGG",
		"AGACCGAGGGAGATAAAGAAATCTAT",
		"CGGGATGATAACCTGCTATGTCACTA",
		"AAACTTTCCAAACGTTTGATATACCC",
	]

	idx = [4, 13, 9, 4, 15]

	sequences_a = ["--" + x + "--" for x in sequences]

	sequences_b = []
	shift = 6
	for i, ix in enumerate(idx):
		seq = sequences[i]
		j = ix - shift
		if j < 0:
			seq = "-"*(-j) + seq
		else:
			seq = seq[j:]
		j = len(sequences[i]) - len(seq)
		if j >= 0:
			seq = seq + "-"*(j)
		else:
			seq = seq[:j]
		sequences_b.append( ("--" + seq)[:25] + "-")

	d = [
		("Rise",  [4 + x for x in [0.97, 0.12, 0.54, 0.33, 0.80, 0.43, 0.50, 0.96, 0.13, 0.91, 0.21]]),
		("Twist", [2 + x for x in [0.70, 0.24, 0.47, 0.07, 0.30, 0.05, 0.85, 0.37, 0.26, 0.49, 0.91]]),
		("Turn",  [0.49, 0.39, 0.86, 0.93, 0.62, 0.95, 0.85, 0.41, 0.10, 0.62, 0.89]),
	]
	dev = [0.24, 0.28, 0.29, 0.27, 0.17, 0.30, 0.29, 0.12, 0.28, 0.26, 0.26]

	fig = plt.figure(figsize=[twocol, onecol])

	plt.subplot(2,2,1)
	ax = plt.gca()
	y_vals = [0, 3, 6, 9, 12]
	for y in y_vals:
		plt.plot(x_vals, helix_a + y, 'steelblue', alpha=0.5)
		plt.plot(x_vals, helix_b + y, 'steelblue', alpha=0.5)


	for i,s in enumerate(sequences_a):
		ax.add_patch(
			mpl.patches.Rectangle(
				(idx[i]*0.9+2.5, y_vals[i]-1), 8, 2, color='darkorange', alpha=0.5)
			)
		ax.text(14, y_vals[i], sequences_a[i], va='center', ha='center', 
			family='monospace', fontsize=11, weight='bold')
	ax.yaxis.set_visible(False)
	ax.xaxis.set_visible(False)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["bottom"].set_visible(False)
	ax.spines["left"].set_visible(False)
	ylim = plt.ylim()
	xlim = plt.xlim()

	plt.subplot(2,2,3)
	ax = plt.gca()
	y_vals = [0, 3, 6, 9, 12]
	for y in y_vals:
		plt.plot(x_vals, helix_a + y, 'steelblue', alpha=0.5)
		plt.plot(x_vals, helix_b + y, 'steelblue', alpha=0.5)

	for i,s in enumerate(sequences_a):
		ax.add_patch(
			mpl.patches.Rectangle(
				(7+2.5, y_vals[i]-1), 8, 2, color='darkorange', alpha=0.5)
			)
		ax.text(14, y_vals[i], sequences_b[i], va='center', ha='center', 
			family='monospace', fontsize=11, weight='bold')
	ax.yaxis.set_visible(False)
	ax.xaxis.set_visible(False)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["bottom"].set_visible(False)
	ax.spines["left"].set_visible(False)

	plt.subplot(2,2,2)
	ax = plt.gca()
	y_vals = [0, 3, 6, 9, 12]

	for i,s in enumerate(pwm_c.split("\n")):
		ax.text(14, y_vals[4 - i], pwm_c.split("\n")[i], va='center', ha='center', 
			family='monospace', fontsize=11, weight='bold')

	plt.ylim(ylim)
	plt.xlim(xlim)

	ax.yaxis.set_visible(False)
	ax.xaxis.set_visible(False)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["bottom"].set_visible(False)
	ax.spines["left"].set_visible(False)

	plt.subplot(2,2,4)

	for i in range(3):
		plt.errorbar(range(len(d[i][1])), d[i][1], yerr=dev, color='steelblue')
		plt.text(-0.5, 0.5 + 2*i, d[i][0], va='center', ha='right',
			family='monospace', fontsize=11, weight='bold')
	plt.axvline(0, color='k')
	plt.plot([0,10],[1.5,1.5], color='k')
	plt.plot([0,10],[3.5,3.5], color='k')
	plt.xlim([-3, 10])
	ax = plt.gca()
	ax.yaxis.set_visible(False)
	ax.xaxis.set_visible(False)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["bottom"].set_visible(False)
	ax.spines["left"].set_visible(False)

	plt.tight_layout()
	return fig

if __name__ == '__main__':
	fig = plotter()
	plt.savefig("figures/figure1.pdf")