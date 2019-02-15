#!/usr/bin/env python

import matplotlib as mpl

mpl.use("Agg")
# mpl.rcParams['pdf.fonttype'] = 42
# mpl.rcParams['font.family'] = 'Myriad Pro'

import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import numpy as np
import string
import sys

from scipy import stats

import overview_fig

np.random.seed(1241)

from Bio import motifs as biomotifs
from strum import strum

######################################################################
# Load the data from all of the outputs                              #
######################################################################

parser = argparse.ArgumentParser(
	description="Generate the figures to accompany the main StruM paper."
	)

parser.add_argument('-a', '--auc', metavar='chip_auc.txt', 
	type=argparse.FileType('rb'), required=True, 
	help="Path to file with area under the ROC curve results for ChIP v Random sequence.")
parser.add_argument('-r', '--correlations', metavar='correlations.txt', 
	type=argparse.FileType('rb'), required=True, 
	help="Path to file with correlations between StruM and PWM scores.")
parser.add_argument('-p', '--positions', metavar='position_comp.txt', 
	type=argparse.FileType('rb'), required=True, 
	help="Path to file with analysis of relative positions between PWM and StruM occurences.")
parser.add_argument('-c', '--coefficients', metavar='coefficients.txt', 
	type=argparse.FileType('rb'), required=True, 
	help="Path to file with coefficients from Logistic Regression model combining PWMs with StruMs.")
parser.add_argument('-s', '--specificities', metavar='specificities.txt', 
    type=argparse.FileType('rb'), required=True, 
    help="Path to file with AUCs of the specificity test.")

args = parser.parse_args()


## AUCs from Peak v. Non-Peak comparison
## Uses 2 sets of background sequences
f = args.auc
pknonpk_labels = ['PWM', 'DWM', 'ML-StruM', 'EM-StruM', 'Combined']
shuff_AUCs, flank_AUCs, pknonpk_genes = [], [], []
shuff_PRCs, flank_PRCs = [], []
while True:
    l1 = f.readline()
    l2 = f.readline()
    if not l2:
        break
    pknonpk_genes.append(l1.split()[0])
    shuff_AUCs.append([float(x) for x in l1.split()[1:6]])
    shuff_PRCs.append([float(x) for x in l1.split()[6:11]])
    flank_AUCs.append([float(x) for x in l2.split()[1:6]])
    flank_PRCs.append([float(x) for x in l2.split()[6:11]])

shuff_AUCs = np.asarray(shuff_AUCs)
shuff_PRCs = np.asarray(shuff_PRCs)
flank_AUCs = np.asarray(flank_AUCs)
flank_PRCs = np.asarray(flank_PRCs)

## Correlations of the scores between each PWMs and EM-StruMs
corr = {}

for line in args.correlations:
	fields = line.split()
	g = fields[0]
	r = float(fields[3])
	corr[g] = r

## Average distance between top StruM and PWM matches
avgdist_tfs, avgdist_ns, avgdist_ks, avgdist_data = [], [], [], []
for line in args.positions:
	fields = line.split()
	avgdist_tfs.append((fields[0]))
	avgdist_ns.append(int(fields[1]))
	avgdist_ks.append(int(fields[2]))
	avgdist_data.append(float(fields[3]))

## Logit model coefficients
f = args.coefficients
coeff_labels = ["PWM", "DWM", "ML-StruM", "EM-StruM"]
coeff_genes = []
shuff_coeffs, flank_coeffs = [], []
while True:
    l1 = f.readline()
    # l2 = f.readline()
    if not l1:
        break
    coeff_genes.append(l1.split()[0])
    shuff_coeffs.append([float(x) for x in l1.split()[1:6]])
    # flank_coeffs.append([float(x) for x in l2.split()[1:6]])

shuff_coeffs = np.asarray(shuff_coeffs)
auc_coeff_idx = {}
for i,g in enumerate(coeff_genes):
	auc_coeff_idx[g] = i
auc_ordr = [auc_coeff_idx[g] for g in pknonpk_genes]
shuff_coeffs = shuff_coeffs[auc_ordr]

## AUCs from TF performance in a TF vs Same-Family and 
## TF vs Other-Family sequences.
###Format: TF auROC-Diff auPRC-Diff auROC-Same auPRC-Same
f = args.specificities
spec_tfs = []
spec_data = []
for line in f:
    fields = line.split()
    spec_tfs.append(fields[0])
    spec_data.append([float(x) for x in fields[1:]])

spec_tfs = np.asarray(spec_tfs)
spec_data = np.asarray(spec_data)

# flank_coeffs = np.asarray(flank_coeffs)

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
	if title is not None:
		plt.title(title,)

	n = data.shape[1]
	positions = range(1, n+1)
	for i in range(n):
		plt.plot(jitter(data[:,i]) + positions[i], data[:,i], '.',
			zorder=1, c=colors[i],
			)

def plot_pvals(p_vals, labels, ax, loc=4):
	cmap_name = "Significance"
	colors = [
			  (255,255,255), 
			  #(239,243,255),
			  (189,215,231),
			  (107,174,214),
			  (33,113,181),]
	colors = [[x/255. for x in y] for y in colors]
	N = len(colors)
	cm = LinearSegmentedColormap.from_list(
        cmap_name, colors, N=N)

	n = p_vals.shape[0]
	for i in range(n):
		for j in range(n):
			if j >= i:
				p_vals[i,j] = np.nan
				continue
			s = get_sig(p_vals[i,j])
			if s == 'ns':
				val = 0
			else:
				val = len(s)
			p_vals[i,j] = val
	axin = inset_axes(ax, width='25%', height='25%', loc=loc, borderpad=2)
	# axin.xaxis.tick_top()
	im = axin.imshow(p_vals, cmap=cm, vmin=-0.5, vmax=3.5)
	plt.yticks(range(1, n), labels[1:])
	plt.xticks(range(n-1), labels[:-1])

	for side in ['top', 'right', 'left', 'bottom']:
		axin.spines[side].set_visible(False)

	cax = inset_axes(
					 axin, width='10%', height='80%', loc=3, 
					 bbox_to_anchor=(0.9, 0, 1, 1), 
					 bbox_transform=axin.transAxes, borderpad=0)
	cbar = plt.colorbar(im, cax=cax, ticks=[0, 1, 2, 3,],)
	cbar.ax.set_yticklabels(['ns', '*', '**', '***'])
	axin.patch.set_visible(False)
	cax.patch.set_visible(False)

def label_plots(fig, plots_high, plots_wide):
	alpha = string.ascii_uppercase
	ref_y = 0.99
	ref_x = 0.01
	del_y = (-1. + 0.05)/plots_high
	del_x = (1. - 0.05)/plots_wide
	for i in range(plots_high):
		for j in range(plots_wide):
			n = i*plots_wide + j
			fig.text(ref_x + del_x*j, ref_y + del_y*i, alpha[n], 
				ha='left', va='top', weight='bold', fontsize=15)

######################################################################
# FUNCTIONS FOR DEALING WITH STATISTICS NICELY
######################################################################

def get_sig(p):
	if 0 <= p < 10.e-10:
		return "***"
	elif 0 <= p < 10.e-5:
		return "**"
	elif 0 <= p < 10.e-2:
		return "*"
	else:
		return "ns"

def print_title(string):
	l = len(string) + 4
	print "\n" + "#"*l
	print "# {} #".format(string.upper())
	print "#"*l

def do_ttest(dataset, labels,):
	n_c = max([len(l) for l in labels])
	header = "{:^{}} v. {:^{}} {:>10} {:>10} {:>10} {:>4}".format("Motif1", n_c, "Motif2", n_c, "Avg diff", "t-value", "p-value", "Sig")
	n = dataset.shape[1]
	l = len(header)
	print "\n" + "="*l
	print header
	print "="*l

	p_vals = np.zeros([dataset.shape[1]]*2)	
	t_vals = np.zeros([dataset.shape[1]]*2)

	for i in range(n):
		for j in range(n):
			t,p = stats.ttest_rel(dataset[:,i], dataset[:,j])
			avgdif = np.average(dataset[:,i] - dataset[:,j])
			p_vals[i,j] = p
			t_vals[i,j] = t
			print "{:^{}} v. {:^{}} {: 10.2e} {: 10.2f} {:10.2e} {:>4}".format(labels[i], n_c, labels[j], n_c,  avgdif, t, p, get_sig(p))
		print "{:^{}}".format("Average {:.2f}".format(np.average(dataset[:,i])), n_c)
		print "-"*l
	return t_vals, p_vals

######################################################################
# Figure 1: Overview figure
######################################################################

fig1 = overview_fig.plotter()
label_plots(fig1, 2, 2)

plt.savefig('figures/figure1.pdf')

######################################################################
# Figure 2: Example motif graphics
######################################################################

# Read example data and construct logos for PWM and StruM for FOXA1
motif = strum.StruM(mode='proteingroove')

head, seq = motif.read_FASTA(open("data/MA0148.1.sites"))
seq = [x.strip(string.ascii_lowercase).replace("X", "N") for x in seq]

motif.train(seq, fasta=False)

motif.PWM *= 100
pwm = biomotifs.read(motif.print_PWM().split('\n'), 'pfm')
pwm.weblogo(
	"output/FOXA1_PWM.png", show_errorbars=False, 
	color_scheme='color_classic', format='png_print')

motif.plot("output/FOXA1_StruM.png")

fig2, (ax_a, ax_b) = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[1, 2]}, 
						   figsize=[onecol, 1.5*onecol])

pwm_img = mpl.image.imread("output/FOXA1_PWM.png")
strum_img = mpl.image.imread("output/FOXA1_StruM.png")
ax_a.imshow(pwm_img)
ax_a.xaxis.set_visible(False)
ax_a.yaxis.set_visible(False)
ax_b.imshow(strum_img)
ax_b.xaxis.set_visible(False)
ax_b.yaxis.set_visible(False)
ax_a.spines["top"].set_visible(False)
ax_a.spines["right"].set_visible(False)
ax_a.spines["bottom"].set_visible(False)
ax_a.spines["left"].set_visible(False)
ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)
ax_b.spines["bottom"].set_visible(False)
ax_b.spines["left"].set_visible(False)

fig2.text(0.01, 0.99, 'A', ha='left', va='top', weight='bold', fontsize=15)
fig2.text(0.01, 0.6, 'B', ha='left', va='top', weight='bold', fontsize=15)

plt.subplots_adjust(left=0.01, bottom=0.05, right=0.99, top=0.95, hspace=0.2)
plt.savefig('figures/figure2.pdf', dpi=600)

######################################################################
# Figure 3: StruMs outperform, but are complementary to, PWMs (AUCs) #
######################################################################

# props = {'color':'white', 'linestyle':'-', 'zorder':9, 'linewidth':3.5}
props2 = {'color':'black', 'linestyle':'-', 'zorder':10, 'linewidth':2}
colors = ["darkorange", "#fb9a99", "seagreen", "steelblue", "mediumpurple"]
markers = ['d', '^', 's', 'o', '*']

fig3 = plt.figure(figsize=[2*twocol, 2*4./5*twocol])
plt.subplot(2,1,1)
label_plots(fig3, 2, 2)

# 3A) auROCs of motifs using shuffled sequence as background
ax_a = plt.subplot(2,2,1)
plt.title("Shuffled background", weight='bold')
plot_box(shuff_AUCs, pknonpk_labels)
plt.ylabel("auROC", weight='bold')

## STATISTICS
print_title("Difference in auROCs, shuffled background")
t_vals, p_vals = do_ttest(shuff_AUCs, pknonpk_labels)
# box_sig_lines(p_vals, ax_a, shuff_AUCs)
plot_pvals(p_vals, [x[0] for x in pknonpk_labels], ax_a)

# 3B) auROCs of motifs using flanking sequence as background
ax_b = plt.subplot(2,2,2)
plt.title("Flanking background", weight='bold')
plot_box(flank_AUCs, pknonpk_labels)
plt.ylabel("auROC", weight='bold')

## STATISTICS
print_title("Difference in auROCs, flanking seq as background")
t_vals, p_vals = do_ttest(flank_AUCs, pknonpk_labels)
# box_sig_lines(p_vals, ax_b, flank_AUCs)
plot_pvals(p_vals, [x[0] for x in pknonpk_labels], ax_b)

# 3C) auPRCs of motifs using shuffled sequence as background
ax_c = plt.subplot(2,2,3)
plot_box(shuff_PRCs, pknonpk_labels)
plt.ylabel("auPRC", weight='bold')

## STATISTICS
print_title("Difference in auPRCs, shuffled background")
t_vals, p_vals = do_ttest(shuff_PRCs, pknonpk_labels)
# box_sig_lines(p_vals, ax_c, shuff_PRCs)
plot_pvals(p_vals, [x[0] for x in pknonpk_labels], ax_c)

# 3D) auPRCs of motifs using flanking sequence as background
ax_d = plt.subplot(2,2,4)
plot_box(flank_PRCs, pknonpk_labels)
plt.ylabel("auPRC", weight='bold')

## STATISTICS
print_title("Difference in auPRCs, flanking seq as background")
t_vals, p_vals = do_ttest(flank_PRCs, pknonpk_labels)
# box_sig_lines(p_vals, ax_d, flank_PRCs)
plot_pvals(p_vals, [x[0] for x in pknonpk_labels], ax_d)

plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.975,
                wspace=0.2, hspace=0.2)

for ax in [ax_a, ax_b, ax_c, ax_d,]:
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)

plt.savefig("figures/figure3.pdf")

######################################################################
# Figure 4: Readout mechanisms revealed by model                     #
######################################################################

fig4 = plt.figure(figsize=[2*twocol, 2*4./5*twocol])
label_plots(fig4, 2, 2)

# 4A) Logit coefficients sorted by improvement over PWM AUC
pwm_auc = shuff_AUCs[:,0]
strum_auc = shuff_AUCs[:,3]

ax_a = plt.subplot(2,2,1)
log_pwm = shuff_AUCs[:,4] - pwm_auc
log_strum = shuff_AUCs[:,4] - strum_auc

print_title("Complementarity of Motifs")
n = np.sum((log_pwm > 0) & (log_strum > 0))
N = len(log_pwm)
print "\nComplementarity for {} experiments".format(n)
print "Fraction: {}".format(float(n)/N)
n2 = np.sum((log_pwm > 0) & (log_strum < 0))
print "\nStruMs are best for {} experiments".format(n2)
print "Fraction: {}".format(float(n2)/N)
print "\nTotal number of experiments: {}".format(N)

plt.plot(log_pwm, log_strum, '.', c='gray')

plt.xlabel("Combined - PWM ($\Delta$ AUC)", weight='bold')
plt.ylabel("Combined - StruM ($\Delta$ AUC)", weight='bold')

# xmin, xmax = plt.xlim()
# ymin, ymax = plt.ylim()
# xmin = ymin = min(xmin, ymin)
# plt.xlim([xmin, xmax])
# plt.ylim([ymin, ymax])
# ax_a.spines["top"].set_visible(False)
# ax_a.spines["right"].set_visible(False)
plt.axhline(0, c='k')
plt.axvline(0, c='k')
ax_a.yaxis.set_label_coords(-0.1,0.5)

for side in ['top', 'right', 'left', 'bottom']:
		ax_a.spines[side].set_visible(False)


# 4B) Specificity of shape- vs base-readout TFs
ax_b = plt.subplot(2,2,2)

# goi_seq = ['STAT', 'GATA']
goi_shp = ['TBP', 'LEF', 'RFX',]
goi_seq = [
	"STAT", "NFAT", "GATA", # Definitely in the dataset
	"BHTH", "homeodomains", # Classes of proteins
		"CTCF", "CTCFL", #"DNMT1", "E4F1", "GATA2", "HINFP", "KLF13", # Zinc fingers
		#"KLF1", "VEZF1", "ZBTB40", "ZFX", "EGR1", "GATA1", "GATAD2A", # Zinc fingers
		#"GATAD2B", "IKZF1", "KLF16", "MTA1", "MYNN", "PRDM10", "REST", # Zinc fingers
		#"RLF", "SP1", "THAP1", "ZBED1", "ZBTB2", "ZBTB33", "ZBTB5", # Zinc fingers
		#"ZBTB7A", "ZFP91", "ZKSCAN1", "ZSCAN29", # Zinc fingers
		# "HMBOX1", "MEIS2", "PKNOX1",  "CUX1", "SIX5", "ZEB2", "ZHX1", # Homeodomains
	"P53", "RAR", # Missing in the dataset
	"NF-kB", "NFRKB", # NF-kB like
	"KORA", # Pseudomondas
	"TALE", # Bacteria
	"TFIIIA", # Arabidopsis
	]
goi_shp = [
	"LEF1", 
	"RFX1", "RFX5", # related to RFX1
	"CMYB", "MYBL2", # related to c-MYB
	"Oct4", "POU5F1",  # Another name for Oct4
	"TBP", "CAP", "RevErb", "SRY", "SOX", # Missing from dataset
	"434", "trp repressor", "MET repressor", "metJ", "IHF", "p22 repressor", "mar A", # Bacteria
	"Scr", "Exd", # Drosophila
	]


text = []

roi_seq = np.zeros(pwm_auc.shape, dtype=bool)
for subg in goi_seq:
	for i, g in enumerate(pknonpk_genes):
		if g.startswith(subg):
			g = g.split('.')[0]
			if len(g) > len(subg) + 1: continue
			roi_seq[i] = True
			text.append((i, g, pwm_auc[i], strum_auc[i], 0))
roi_shp = np.zeros(pwm_auc.shape, dtype=bool)
for subg in goi_shp:
	for i, g in enumerate(pknonpk_genes):
		if g.startswith(subg):
			g = g.split('.')[0]
			if len(g) > len(subg) + 1: continue
			roi_shp[i] = True
			text.append((i, g, pwm_auc[i], strum_auc[i], 1))

plt.plot(pwm_auc, strum_auc, '.', c='gray')
plt.plot(pwm_auc[roi_seq], strum_auc[roi_seq], markers[0], 
		 c=colors[0], ms=6, label='Base Readout')
plt.plot(pwm_auc[roi_shp], strum_auc[roi_shp], markers[3], 
		 c=colors[3], ms=6, label='Shape Readout')
plt.plot(plt.xlim(), plt.ylim(), '--', c='gray')

for i,g,x,y,j in text:
	plt.text(x, y, g, dict(weight='bold'))

plt.legend(loc='upper left')
plt.xlabel("PWM auROC", weight='bold')
plt.ylabel("EM-StruM auROC", weight='bold')

extra_colors = ['gold', 'skyblue']

# 4C)  Logit coefficients sorted by improvement over PWM AUC
ax_c = plt.subplot(2,2,3)

plt.ylabel("logit coefficients", weight='bold')
x = shuff_AUCs[:,4] - shuff_AUCs[:,0]
idx = np.argsort(x)
for i in [0,3]:#range(4):
	# plt.plot(x, shuff_coeffs[:,i], markers[i], c=colors[i], label=coeff_labels[i])
	plt.plot(shuff_coeffs[:,i][idx], markers[i], c=colors[i], label=coeff_labels[i])
for i,g,xval,yval,j in text:
	xval = np.where(idx==i)[0]
	y1 = shuff_coeffs[i,0]
	y2 = shuff_coeffs[i,3]
	yt = (y1+y2)/2.
	plt.plot([xval,xval], [y1,y2], '-', color=extra_colors[j], linewidth=2.5)
	plt.text(xval, yt, g, dict(weight='bold'))
plt.legend(ncol=4, bbox_to_anchor=[0.0, 0.93], loc='upper left')
ax_c2 = ax_c.twinx()
ax_c2.tick_params('y', colors='r')
plt.plot(x[idx], 'r-', label="Combined - PWM ($\Delta$ AUC)")
plt.legend(loc='upper left', bbox_to_anchor=[0.0, 1])
ax_c.set_xticks([])
ax_c.set_xlabel("Ranked by $\Delta$ AUC", weight='bold')

ax_c.yaxis.set_label_coords(-0.1,0.5)
ax_c2.yaxis.set_label_coords(1.15,0.5)
# plt.xlabel("Combined - PWM ($\Delta$ AUC)", weight='bold')

# 4D) Logit coefficients sorted by improvement over StruM AUC
ax_d = plt.subplot(2,2,4)
# plt.ylabel("logit coefficients", weight='bold')
x = shuff_AUCs[:,4] - shuff_AUCs[:,3]
idx = np.argsort(x)
for i in [0,3]:#range(4):
	# plt.plot(x, shuff_coeffs[:,i], markers[i], c=colors[i], label=coeff_labels[i])
	plt.plot(shuff_coeffs[:,i][idx], markers[i], c=colors[i], label=coeff_labels[i])
for i,g,xval,yval,j in text:
	xval = np.where(idx==i)[0]
	y1 = shuff_coeffs[i,0]
	y2 = shuff_coeffs[i,3]
	yt = (y1+y2)/2.
	plt.plot([xval,xval], [y1,y2], '-', color=extra_colors[j], linewidth=2.5)
	plt.text(xval, yt, g, dict(weight='bold'))
ax_d2 = ax_d.twinx()
# ax_d2.set_ylabel("Combined - StruM ($\Delta$ AUC)", color='r', weight='bold', rotation=270)
ax_d2.tick_params('y', colors='r')
plt.plot(x[idx], 'r-', label="Combined - StruM ($\Delta$ AUC)")
plt.legend(loc='upper left', bbox_to_anchor=[0.0, 1])
ax_d.set_xticks([])
ax_d.set_xlabel("Ranked by $\Delta$ AUC", weight='bold')
ax_d.yaxis.set_label_coords(-0.1,0.5)
ax_d2.yaxis.set_label_coords(1.15,0.5)

# plt.legend(ncol=4, bbox_to_anchor=[0.5, 1], loc='upper center')
# plt.xlabel("Combined - StruM ($\Delta$ AUC)", weight='bold')

plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.975,
                wspace=0.2, hspace=0.2)

for ax in [ax_a, ax_b, ax_c, ax_d,]:
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)

plt.savefig("figures/figure4.pdf")


######################################################################
# Figure 5: Distribution of significant StruM matches vs PWM matches #
######################################################################

fig5 = plt.figure(figsize=[2*onecol, 2.5*twocol])
label_plots(fig5, 5, 1)

ax_a = plt.subplot(5,1,1)
plt.hist(avgdist_data, bins=np.linspace(0,60,31))

# Correlated
ax_b = plt.subplot(5,1,2)
cor_img = mpl.image.imread("output/" + "TAL1.ENCFF519DOC" # "eGFP-GTF2E2.ENCFF394WVZ" # "FOSL1.ENCFF961OPH"
							+ "_fimo_v_strum_matches2.png")
ax_b.set_title("TAL1", loc='left')
ax_b.imshow(cor_img)
ax_b.xaxis.set_visible(False)
ax_b.yaxis.set_visible(False)

# Flanking
ax_c = plt.subplot(5,1,3)
flank_img = mpl.image.imread("output/" + "RFX1.ENCFF934JXG" # "eGFP-ZNF512.ENCFF617CTX" # ZBTB33.ENCFF681IOP
								+ "_fimo_v_strum_matches2.png")
ax_c.set_title("RFX1", loc='left')
ax_c.imshow(flank_img)
ax_c.xaxis.set_visible(False)
ax_c.yaxis.set_visible(False)


# Anti Correlated
ax_d = plt.subplot(5,1,4)
anticor_img = mpl.image.imread("output/" + "MEF2A.ENCFF883WDT" # "RLF.ENCFF569QYK" # "ATF4.ENCFF491DNM" # ATF2.ENCFF525YRJ
							+ "_fimo_v_strum_matches2.png")
ax_d.set_title("RLF", loc='left')
ax_d.imshow(anticor_img)
ax_d.xaxis.set_visible(False)
ax_d.yaxis.set_visible(False)

# Not Correlated
ax_e = plt.subplot(5,1,5)
rand_img = mpl.image.imread("output/" + "SNIP1.ENCFF772GVZ" # "ATF4.ENCFF491DNM" # ATF2.ENCFF525YRJ
							+ "_fimo_v_strum_matches2.png")
ax_e.set_title("SNIP1", loc='left')
ax_e.imshow(rand_img)
# ax_e.xaxis.set_visible(False)
ax_e.yaxis.set_visible(False)
xmin, xmax = ax_e.get_xlim()
xticks = range(int(xmin), int(xmax+1), int((xmax-xmin)/4))
xticks = np.linspace(xmin, xmax, 5)
ax_e.set_xticks(xticks)
ax_e.set_xticklabels(['-100', '0\nPWM Position (bp)', '+/-100', '0\nStruM Position (bp)', '+100'])

plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.975,
                wspace=0.2, hspace=0.2)
plt.savefig('figures/figure5.pdf')

################################################
# Figure 6: Correlation between PWM and StruMs #
################################################

R = []
for g in pknonpk_genes:
	if g in corr:
		R.append(corr[g])
left = 0.2
bottom = 0.2

fig5 = plt.figure(figsize=[onecol, onecol])
label_plots(fig5, 1, 1)
x = R
y = log_pwm
ax_a = fig5.add_axes((left, bottom, 0.8-left, 0.8-bottom))
ax_pos = ax_a.get_position().bounds
ax_top = fig5.add_axes((left, 0.8, ax_pos[2], 0.2), sharex=ax_a)
ax_right = fig5.add_axes((0.8, bottom, 0.2, ax_pos[3]), sharey=ax_a)
ax_a.scatter(x, y, s=5, c='k')
ax_a.set_xlabel("Correlation between scores", weight='bold')
ax_a.set_ylabel("$\\Delta$ AUC (logit-PWM)", weight='bold')
ax_top.hist(x, bins=30, facecolor='grey')
ax_right.hist(y, bins=30, facecolor='grey', orientation='horizontal')
ax_top.axis('off')
ax_right.axis('off')
z = np.polyfit(x, y, 1)
p = np.poly1d(z)
minx = min(x)
maxx = max(x)
ax_a.plot([minx,maxx],p([minx,maxx]), c=colors[0], lw=2)
fig5.savefig('figures/figure6.pdf')

print """
############################################
# CORRELATION BETWEEN PWM AND STRUM SCORES #
############################################
"""
print "Eqn of the line:   y = {:0.3f} x + {:0.3f}".format(z[0], z[1])
print "Correlation: {:0.3f}".format(np.corrcoef(x,y)[0,1])

#####################################
# Figure 7: Specificities of StruMs #
#####################################

fig6 = plt.figure(figsize([twocol, onecol]))
label_plots(fig6, 1, 2)

# auROC
ax_a =  fig6.add_subplot(1, 2, 1)
x1 = spec_data[:,0]
y1 = spec_data[:,2]
minx = min(np.min(x1), np.min(y1))
maxx = max(np.max(x1), np.max(y1))
ax_a.plot(x1, y1, '.', c=colors[3])
ax_a.plot([minx, maxx], [minx, maxx], '--', c='gray')
ax_a.set_xlim([minx, maxx])
ax_a.set_ylim([minx, maxx])
ax_a.set_xlabel("AUC (TF Peak vs Other TF Families)", weight='bold')
ax_a.set_ylabel("AUC (Same TF family vs Other TF Families)", weight='bold')

# auPRC
ax_b =  fig6.add_subplot(1, 2, 2)
x2 = spec_data[:,1]
y2 = spec_data[:,3]
minx = min(np.min(x2), np.min(y2))
maxx = max(np.max(x2), np.max(y2))
ax_b.plot(x2, y2, '.', c=colors[3])
ax_b.plot([minx, maxx], [minx, maxx], '--', c='gray')
ax_b.set_xlim([minx, maxx])
ax_b.set_ylim([minx, maxx])
ax_b.set_xlabel("AUC (TF Peak vs Other TF Families)", weight='bold')
ax_b.set_ylabel("AUC (Same TF family vs Other TF Families)", weight='bold')