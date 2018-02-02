#!/usr/bin/env python

import matplotlib as mpl
# mpl.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.family'] = 'Myriad Pro'
import matplotlib.pyplot as plt
import numpy as np
import sys
# plt.style.use("dark_background")

f = open(sys.argv[1])

simple_shuff, flanking, genes = [], [], []

while True:
    l1 = f.readline()
    l2 = f.readline()
    if not l2:
        break
    genes.append(l1.split()[0])
    simple_shuff.append([float(x) for x in l1.split()[1:6]])
    flanking.append([float(x) for x in l2.split()[1:6]])

simple_shuff = np.asarray(simple_shuff)
flanking = np.asarray(flanking)

def jitter(row, nbins=30, binsize=None, vmin=None, vmax=None):
	if vmin is None:
		minimum = np.min(row)
	else:
		minimum = vmin
	if vmax is None:
		maximum = np.max(row)
	else:
		maximum = vmax
	if binsize is None:
		binsize = (maximum - minimum)/float(nbins)

	bins = np.asarray(row // binsize, dtype=int)
	bins -= np.min(bins)
	bincounts = np.bincount(bins)
	count = np.max(bincounts)
	if count % 2 == 0:
		count += 1
	width = min(0.8, 0.05*count)

	pool = sorted(np.linspace(-width/2, width/2, count), key=lambda x:abs(x))
	xs = np.zeros(row.shape)
	for i in range(bincounts.shape[0]):
		xs[bins == i] = np.random.permutation(pool[:bincounts[i]])
	return xs

plt.figure(figsize=[12,5])

props = {'color':'white', 'linestyle':'-', 'zorder':9, 'linewidth':3.5}
props2 = {'color':'black', 'linestyle':'-', 'zorder':10, 'linewidth':2}
colors = ["darkorange", "#fb9a99", "seagreen", "steelblue", "mediumpurple"]
labels = ['PWM', 'DWM', 'ML-StruM', 'EM-StruM', 'Combined']

def plot_box(data, title, labels,):
	plt.boxplot(data, notch=False, boxprops=props, whiskerprops=props, medianprops=props, sym='')
	plt.boxplot(data, notch=False, boxprops=props2, whiskerprops=props2, medianprops=props2, sym='')
	plt.xticks(range(1, len(labels)+1), labels,)
	plt.ylim([0,1])
	plt.ylabel("AUC",)
	plt.title(title,)

	n = data.shape[1]
	positions = range(1, n+1)
	for i in range(n):
		plt.plot(jitter(data[:,i]) + positions[i], data[:,i], '.',
			zorder=1, c=colors[i], #alpha=0.5
			)


plt.subplot(121)
plot_box(simple_shuff, "Simple shuffle", labels)

plt.subplot(122)
plot_box(flanking, "Flanking Sequences", labels)


plt.tight_layout()

plt.savefig('figures/performance.pdf')
plt.close()

plt.figure(figsize=[10,5])
plt.subplot(121)
plot_box(simple_shuff[:,:-1], "Simple shuffle", labels[:-1])
plt.subplot(122)
plot_box(flanking[:,:-1], "Flanking Sequences", labels[:-1])
plt.tight_layout()
plt.savefig('figures/performance_sub.pdf')
plt.close()

############################################################
# PLOT 2
############################################################
pwm = simple_shuff[:,0]
dwm = simple_shuff[:,1]
ml_strum = simple_shuff[:,2]
em_strum = simple_shuff[:,3]
logit = simple_shuff[:,4]

goi_seq = ['STAT', 'GATA']
goi_shp = ['TBP', 'LEF', 'RFX',]
text = []
roi_seq = np.zeros(pwm.shape, dtype=bool)
for subg in goi_seq:
	for i, g in enumerate(genes):
		if g.startswith(subg):
			roi_seq[i] = True
			text.append((g, pwm[i], em_strum[i]))
roi_shp = np.zeros(pwm.shape, dtype=bool)
for subg in goi_shp:
	for i, g in enumerate(genes):
		if g.startswith(subg):
			roi_shp[i] = True
			text.append((g, pwm[i], em_strum[i]))

plt.figure(figsize=[6,6])
plt.plot(pwm, em_strum, '.', c='gray')
plt.plot(pwm[roi_seq], em_strum[roi_seq], 'd', c=colors[3], ms=12, label='Base Readout')
plt.plot(pwm[roi_shp], em_strum[roi_shp], 's', c=colors[0], ms=12, label='Shape Readout')
plt.plot([np.min(pwm), np.max(pwm)], [np.min(em_strum), np.max(em_strum)], '--', c='gray')

for g,x,y in text:
	plt.text(x, y, g, dict(weight='bold'))

plt.legend()

plt.xlabel("PWM auROC")
plt.ylabel("EM-StruM auROC")
plt.tight_layout()
plt.savefig("figures/pwm_v_strum_auc.pdf")
plt.close()

############################################################
# PLOT 3
############################################################

plt.figure(figsize=[6,6])
plt.axvline(0)
plt.axhline(0)

delta = simple_shuff[:,[0,3]]
delta[:,0] = simple_shuff[:,4] - delta[:,0]
delta[:,1] = simple_shuff[:,4] - delta[:,1]

pos_delta = delta[np.min(delta, axis=1) > 0]
neg_delta = delta[np.min(delta, axis=1) <= 0]

plt.plot(delta[:,0], delta[:,1], '.', c='gray')

plt.xlabel("Combined - PWM\n($\Delta$ AUC)")
plt.ylabel("Combined - StruM\n($\Delta$ AUC)")

xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
xmin = ymin = min(-0.1, xmin, ymin)
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])

ax = plt.gca()

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)

plt.tight_layout()
plt.savefig("figures/logit_improvement.pdf")
plt.close()

############################################################
# T TESTS
############################################################

from scipy.stats import ttest_rel

def ttest(X):
	header = "| {:^8s} | {:<13s} | ".format("Model", "Avg AUC (SD)")
	header2 = " | ".join(["{:<8s}".format(l) for l in labels])
	divider = "|:--------:|:--------------" + "|:---------"*len(labels) + "|"
	print header + header2 + " |"
	print divider
	for i in range(X.shape[1]):
		avg = np.average(X[:,i])
		std = np.std(X[:,i])
		pees = []
		for j in range(X.shape[1]):
			t,p = ttest_rel(X[:,i], X[:,j])
			pees.append(p)

		nl = "| {:^8s} | {:5.3f} ({:5.3f}) | ".format(labels[i], avg, std)
		nl2 = " | ".join(["{:<8.2e}".format(x) for x in pees])
		print nl + nl2 + " |"

print "\nn = {}\n".format(simple_shuff.shape[0])

print "\n##Shuffled Background\n"
ttest(simple_shuff)

print "\n##Flanking Sequences\n"
ttest(flanking)


############################################################
# CORRELATION V AUC
############################################################

f = open(sys.argv[2])

corr = {}

for line in f:
	fields = line.split()
	g = fields[0]
	r = float(fields[3])
	corr[g] = r

x,y,c,s = [],[],[],[]
for i,g in enumerate(genes):
	if g in corr:
		y.append(abs(corr[g]))
		x.append(simple_shuff[i,3] - simple_shuff[i,0])
		for gsub in goi_seq:
			if g.startswith(gsub):
				c.append('steelblue')
				s.append(12)
				break
		else:
			for gsub in goi_shp:
				if g.startswith(gsub):
					c.append('darkorange')
					s.append(12)
					break
			else:
				c.append('gray')
				s.append(8)

plt.figure(figsize=[6,6])
plt.scatter(x,y,c=c,s=s)
plt.xlabel("EM-StruM - PWM\n($\Delta$ AUC)")
plt.ylabel("Correlation $|r|$")
plt.savefig("figures/corr_perf.pdf")
plt.close()

############################################################
# TRY TO COMPARE TO DNASE
############################################################
"""
cd ~/Desktop/StruM/performance/dnase
ssh pdeford1@comp1.bx.bio.jhu.edu

cd strum_proof/DNase

ls output/*_AUCs.txt | while read line;
do
	fname=${line##*/}
	name=${fname%%.*}
	tf=${name%%_*}
	echo -n $tf$'\t'
	awk '{printf "%s\t",substr($NF,2,4)}' 'output/'$tf'_AUCs.txt'
	echo
done > output/dnase_consolidated.txt

logout

scp pdeford1@comp1.bx.bio.jhu.edu:strum_proof/DNase/output/dnase_consolidated.txt ./output
"""

x,y,c = [],[],[]

f = open(sys.argv[3])
for line in f:
	fields = line.split()
	g = fields[0]
	auc = float(fields[-1])
	for i,g2 in enumerate(genes):
		if g2.startswith(g):
			x.append(simple_shuff[i,3])
			y.append(auc)
			c.append(simple_shuff[i,3]-simple_shuff[i,0])

plt.figure(figsize=[6,5])
plt.scatter(x, y, cmap='PiYG', c=c, )#edgecolor='white')#'BrBG'

xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
xmin = ymin = min(xmin, ymin)
xmax = ymax = max(xmax, ymax)
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])

plt.ylabel("DNase StruM AUC\n(Bound $v$ Unbound)")
plt.xlabel("EM-StruM AUC\n(Peak $v$ Non)")

plt.colorbar()
plt.tight_layout()

plt.savefig("figures/dnase_v_em.pdf")
plt.close()


############################################################
# Summarize average distance between top `n` StruM matches
# to the nearest of the top `n` PWM matches

f = open(sys.argv[4])

tfs = []
ns = []
ks = []
dists = []


for line in f:
	fields = line.split()
	tf = (fields[0])
	n = int(fields[1])
	k = int(fields[2])
	dist = float(fields[3])

	tfs.append(tf)
	ns.append(n)
	ks.append(k)
	dists.append(dist)

plt.figure(figsize=[6,4])
plt.hist(dists, bins=30)
plt.xlabel("Average Distance to Nearest FIMO Match (bp)")
ax = plt.gca()
ax.tick_params(labelleft='off')
plt.tight_layout()
plt.savefig("figures/avg_distance_dist.pdf")
