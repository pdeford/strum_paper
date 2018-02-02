#!/usr/bin/env python

import matplotlib as mpl

mpl.use("Agg")
# mpl.rcParams['pdf.fonttype'] = 42
# mpl.rcParams['font.family'] = 'Myriad Pro'

import matplotlib.pyplot as plt
import numpy as np
import string
import sys

np.random.seed(1241)

######################################################################
# Load the data from all of the outputs                              #
######################################################################

## AUCs from Peak v. Non-Peak comparison
## Uses 2 sets of background sequences
f = open(sys.argv[1])
pknonpk_labels = ['PWM', 'DWM', 'ML-StruM', 'EM-StruM', 'Combined']
shuff_AUCs, flank_AUCs, pknonpk_genes = [], [], []
while True:
    l1 = f.readline()
    l2 = f.readline()
    if not l2:
        break
    pknonpk_genes.append(l1.split()[0])
    shuff_AUCs.append([float(x) for x in l1.split()[1:6]])
    flank_AUCs.append([float(x) for x in l2.split()[1:6]])

shuff_AUCs = np.asarray(shuff_AUCs)
flank_AUCs = np.asarray(flank_AUCs)

## Correlations of the scores between each PWMs and EM-StruMs
f = open(sys.argv[2])

corr = {}

for line in f:
	fields = line.split()
	g = fields[0]
	r = float(fields[3])
	corr[g] = r

## AUCs for Performance of Cell type specific predictions (w/ DNase)
f = open(sys.argv[3])
celspec_labels = np.asarray(
	["PWM", 
	 "PWM + Millipede", 
	 "StruM", 
	 "StruM + Millipede",
	 "DNase-StruM"]
)
celspec_data = []
celspec_genes = []
for line in f:
	fields = line.split()
	celspec_genes.append(fields[0])
	celspec_data.append([float(x) for x in fields[1:]])

celspec_data = np.asarray(celspec_data)
sorter = [0,2,1,3,4]
celspec_data = celspec_data[:, sorter]
celspec_labels = celspec_labels[sorter]

## Average distance between top StruM and PWM matches
f = open(sys.argv[4])
avgdist_tfs, avgdist_ns, avgdist_ks, avgdist_data = [], [], [], []
for line in f:
	fields = line.split()
	avgdist_tfs.append((fields[0]))
	avgdist_ns.append(int(fields[1]))
	avgdist_ks.append(int(fields[2]))
	avgdist_data.append(float(fields[3]))

## Logit model coefficients
f = open(sys.argv[5])
coeff_labels = ["PWM", "DWM", "ML-StruM", "EM-StruM"]
coeff_genes = []
shuff_coeffs, flank_coeffs = [], []
while True:
    l1 = f.readline()
    l2 = f.readline()
    if not l2:
        break
    coeff_genes.append(l1.split()[0])
    shuff_coeffs.append([float(x) for x in l1.split()[1:6]])
    flank_coeffs.append([float(x) for x in l2.split()[1:6]])

shuff_coeffs = np.asarray(shuff_coeffs)
flank_coeffs = np.asarray(flank_coeffs)

######################################################################
# Functions to make plotting figures easier
######################################################################

def cm2inch(value):
    return value/2.54

onecol = cm2inch(84/10.)
twocol = cm2inch(178/10.)
maxheight = cm2inch(230/10.)

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

def plot_box(data, labels, title=None):
	plt.boxplot(data, notch=False, boxprops=props2, whiskerprops=props2, medianprops=props2, sym='')
	plt.xticks(range(1, len(labels)+1), labels,)
	plt.ylim([0,1])
	plt.ylabel("AUC", weight='bold')
	if title is not None:
		plt.title(title,)

	n = data.shape[1]
	positions = range(1, n+1)
	for i in range(n):
		plt.plot(jitter(data[:,i]) + positions[i], data[:,i], '.',
			zorder=1, c=colors[i],
			)

def label_plots(fig, plots_high, plots_wide):
	alpha = string.ascii_uppercase
	ref_y = 0.99
	ref_x = 0.01
	del_y = (-1. + 0.05)/plots_high
	del_x = 1./plots_wide
	for i in range(plots_high):
		for j in range(plots_wide):
			n = i*plots_wide + j
			fig.text(ref_x + del_x*j, ref_y + del_y*i, alpha[n], 
				ha='left', va='top', weight='bold', fontsize=15)

######################################################################
# Figure 1: Overview figure
######################################################################

plt.figure(figsize=[onecol, 2*onecol])
ax_a = plt.subplot(4,1,1)

ax_b = plt.subplot(4,1,2)

ax_c = plt.subplot(4,1,3)

ax_d = plt.subplot(4,1,4)

plt.savefig('figures/figure1.pdf')

######################################################################
# Figure 2: Example motif graphics
######################################################################

plt.figure(figsize=[onecol, 2*onecol])
ax_a = plt.subplot(2,1,1)

ax_b = plt.subplot(2,1,2)

plt.savefig('figures/figure2.pdf')

######################################################################
# Figure 3: StruMs outperform, but are complementary to, PWMs        #
######################################################################

# props = {'color':'white', 'linestyle':'-', 'zorder':9, 'linewidth':3.5}
props2 = {'color':'black', 'linestyle':'-', 'zorder':10, 'linewidth':2}
colors = ["darkorange", "#fb9a99", "seagreen", "steelblue", "mediumpurple"]
markers = ['d', '^', 's', 'o', '*']

fig3 = plt.figure(figsize=[2*twocol, 2*2/3.*twocol])
label_plots(fig3, 2, 3)

# 3A) AUCs of motifs using shuffled sequence as background
ax_a = plt.subplot(2,3,1)
plot_box(shuff_AUCs, pknonpk_labels)


# 3B) AUCs of motifs using flanking sequence as background
ax_b = plt.subplot(2,3,2)
plot_box(flank_AUCs, pknonpk_labels)


# 3C) Logit coefficients sorted by improvement over PWM AUC
ax_c = plt.subplot(2,3,3)
x = shuff_AUCs[:,-1] - shuff_AUCs[:,0]
for i in [0,3]:#range(4):
	plt.plot(x, shuff_coeffs[:,i], markers[i], c=colors[i], label=coeff_labels[i])
plt.xlabel("Combined - PWM ($\Delta$ AUC)", weight='bold')

# 3D) Specificity of shape- vs base-readout TFs
ax_d = plt.subplot(2,3,4)

goi_seq = ['STAT', 'GATA']
goi_shp = ['TBP', 'LEF', 'RFX',]
text = []

pwm_auc = shuff_AUCs[:,0]
strum_auc = shuff_AUCs[:,3]

roi_seq = np.zeros(pwm_auc.shape, dtype=bool)
for subg in goi_seq:
	for i, g in enumerate(pknonpk_genes):
		if g.startswith(subg):
			roi_seq[i] = True
			text.append((g, pwm_auc[i], strum_auc[i]))
roi_shp = np.zeros(pwm_auc.shape, dtype=bool)
for subg in goi_shp:
	for i, g in enumerate(pknonpk_genes):
		if g.startswith(subg):
			roi_shp[i] = True
			text.append((g, pwm_auc[i], strum_auc[i]))

plt.plot(pwm_auc, strum_auc, '.', c='gray')
plt.plot(pwm_auc[roi_seq], strum_auc[roi_seq], markers[0], 
		 c=colors[0], ms=6, label='Base Readout')
plt.plot(pwm_auc[roi_shp], strum_auc[roi_shp], markers[3], 
		 c=colors[3], ms=6, label='Shape Readout')
plt.plot(plt.xlim(), plt.ylim(), '--', c='gray')

for g,x,y in text:
	plt.text(x, y, g, dict(weight='bold'))

plt.legend(loc='upper left')
plt.xlabel("PWM auROC", weight='bold')
plt.ylabel("EM-StruM auROC", weight='bold')

# 3E) Complementarity of motifs in LOGIT model
ax_e = plt.subplot(2,3,5)

log_pwm = shuff_AUCs[:,4] - pwm_auc
log_strum = shuff_AUCs[:,4] - strum_auc

plt.plot(log_pwm, log_strum, '.', c='gray')

plt.xlabel("Combined - PWM ($\Delta$ AUC)", weight='bold')
plt.ylabel("Combined - StruM ($\Delta$ AUC)", weight='bold')

xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
xmin = ymin = min(xmin, ymin)
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])


ax_e.spines["top"].set_visible(False)
ax_e.spines["right"].set_visible(False)

# 3F) Logit coefficients sorted by improvement over PWM AUC
ax_f = plt.subplot(2,3,6)
x = shuff_AUCs[:,-1] - shuff_AUCs[:,3]
for i in [0,3]:#range(4):
	plt.plot(x, shuff_coeffs[:,i], markers[i], c=colors[i], label=coeff_labels[i])
plt.legend(ncol=4, bbox_to_anchor=[0.5, 1], loc='upper center')
plt.xlabel("Combined - StruM ($\Delta$ AUC)", weight='bold')

plt.subplots_adjust(left=0.05, bottom=0.1, right=0.99, top=0.975,
                wspace=0.2, hspace=0.2)

for ax in [ax_a, ax_b, ax_c, ax_d, ax_e, ax_f ]:
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)

plt.savefig("figures/figure3.pdf")

######################################################################
# Figure 4: Distribution of significant StruM matches vs PWM matches #
######################################################################

plt.figure(figsize=[2*twocol, 2*2*twocol])

plt.savefig('figures/figure4.pdf')


######################################################################
# Figure 5: Cell-type Specific Prediction Performance                #
######################################################################
colors = ["darkorange", "steelblue", "darkorange", "steelblue", "seagreen"]

plt.figure(figsize=[3*onecol, 2*onecol])

plot_box(celspec_data, celspec_labels)

plt.savefig('figures/figure5.pdf')