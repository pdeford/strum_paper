#!/usr/bin/env python

import matplotlib as mpl
mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

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

# for thing in sequences_a + sequences_b:
# 	print thing

# print " "*8 + "-"*9 + " "*9 + "|"

pwm_c = """   1   2   3   4   5   6   7   8   9
A 0.2 0.7 0.0 1.0 0.0 1.0 0.9 0.2 0.4
C 0.4 0.0 0.0 0.0 0.0 0.0 0.0 0.2 0.2
G 0.3 0.0 1.0 0.0 0.0 0.0 0.1 0.6 0.3
T 0.1 0.3 0.0 0.0 1.0 0.0 0.0 0.0 0.1"""

d = [
	("Rise",  [0.97, 0.12, 0.54, 0.33, 0.80, 0.43, 0.50, 0.96, 0.13, 0.91, 0.21]),
	("Twist", [0.70, 0.24, 0.47, 0.07, 0.30, 0.05, 0.85, 0.37, 0.26, 0.49, 0.91]),
	("Turn",  [0.49, 0.39, 0.86, 0.93, 0.62, 0.95, 0.85, 0.41, 0.10, 0.62, 0.89]),
]
dev = [0.24, 0.28, 0.29, 0.27, 0.17, 0.30, 0.29, 0.12, 0.28, 0.26, 0.26]

plt.figure()
for i in range(3):
	plt.subplot(3,1,i+1)
	plt.ylabel(d[i][0])
	# plt.plot()
	plt.errorbar(range(len(d[i][1])), d[i][1], yerr=dev)
	plt.xticks([])
	plt.yticks([])
plt.xlabel("Position")

plt.savefig("test.png")