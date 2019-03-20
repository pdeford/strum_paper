#!/usr/bin/env python

import matplotlib as mpl

mpl.use("Agg")
# mpl.rcParams['pdf.fonttype'] = 42
# mpl.rcParams['font.family'] = 'Myriad Pro'

import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patheffects as pe

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
pknonpk_labels = ['PWM', 'DWM', 'ML\nStruM', 'EM\nStruM', 'Combined']
pknonpk_labels2 = ['PWM', 'DWM', 'ML-StruM', 'EM-StruM', 'Combined']
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
coeff_labels = ["PWM", "EM-StruM"]
coeff_genes = []
shuff_coeffs, flank_coeffs = [], []
while True:
    l1 = f.readline()
    # l2 = f.readline()
    if not l1:
        break
    coeff_genes.append(l1.split()[0])
    shuff_coeffs.append([float(x) for x in l1.split()[1:]])

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

plt.rcParams.update({'font.size': 10, 'font.family':'sans-serif', 'font.sans-serif':['Arial', 'Helvetica']})
# props = {'color':'white', 'linestyle':'-', 'zorder':9, 'linewidth':3.5}
props2 = {'color':'black', 'linestyle':'-', 'zorder':10, 'linewidth':2}
colors = ["darkorange", "#fb9a99", "seagreen", "steelblue", "mediumpurple"]
markers = ['d', '^', 's', 'o', '*']

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
            zorder=1, c=colors[i], markersize=3,
            )

def plot_pvals(p_vals, labels, ax, loc=4, cbar=True):
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
    axin = inset_axes(ax, width='100%', height='100%', 
                      loc=loc, borderpad=0,
                      bbox_to_anchor=(0.8, 0.15, 0.25, 0.25), 
                      bbox_transform=ax.transAxes, )
    # axin.xaxis.tick_top()
    im = axin.imshow(p_vals, cmap=cm, vmin=-0.5, vmax=3.5)
    plt.yticks(range(1, n), labels[1:], ha='center', size='small')
    plt.xticks(range(n-1), labels[:-1], size='small')

    for side in ['top', 'right', 'left', 'bottom']:
        axin.spines[side].set_visible(False)

    if cbar:
        cax = inset_axes(
                         axin, width='10%', height='80%', loc=3, 
                         bbox_to_anchor=(0.9, 0, 1, 1), 
                         bbox_transform=axin.transAxes, borderpad=0)
        cbar = plt.colorbar(im, cax=cax, ticks=[0, 1, 2, 3,],)
        cbar.ax.set_yticklabels(['ns', '*', '**', '***'])
        cax.patch.set_visible(False)
    axin.patch.set_visible(False)

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
                ha='left', va='top', weight='bold', fontsize=10)

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

print_title('PWM FOR FOXA1 (MA0148.1)')

motif.PWM *= 100
pwm = biomotifs.read(motif.print_PWM().split('\n'), 'pfm')
##### UNCOMMENT THIS FROM HERE #####
pwm.weblogo(
    "output/FOXA1_PWM.png", show_errorbars=False, 
    color_scheme='color_classic', format='png_print')

motif.plot("output/FOXA1_StruM.png")
##### TO HERE #####

# fig2, (ax_a, ax_b) = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[1, 4]}, 
#                            figsize=[onecol, 1.2*onecol])

fig2 = plt.figure(figsize=(onecol, 1.3*onecol))
ax_a = fig2.add_axes((0.11, 0.67, 0.95, 0.325))
ax_b = fig2.add_axes((0.0,   0.0, 1.0, 0.7))

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

fig2.text(0.01, 0.99, 'A', ha='left', va='top', weight='bold', fontsize=10)
fig2.text(0.01, 0.7, 'B', ha='left', va='top', weight='bold', fontsize=10)

plt.subplots_adjust(left=0.01, bottom=0.05, right=0.99, top=0.95, hspace=0.2)
plt.savefig('figures/figure2.pdf', dpi=600)

#####################################
# Figure 3: Specificities of StruMs #
#####################################

fig3 = plt.figure(figsize=(onecol, 0.8*onecol))
label_plots(fig3, 1, 2)

# auROC
ax_a =  fig3.add_subplot(1, 2, 1)
x1 = spec_data[:,0]
y1 = spec_data[:,2]
plot_box(data=np.vstack([x1,y1]).T, labels=["Cross\nFamily", "Control"], title="auROC")
# for tick in ax_a.get_xticklabels():
#     tick.set_rotation(20)

# auPRC
ax_b =  fig3.add_subplot(1, 2, 2)
x2 = spec_data[:,1]
y2 = spec_data[:,3]
plot_box(data=np.vstack([x2,y2]).T, labels=["Cross\nFamily", "Control"], title="auPRC")
# for tick in ax_b.get_xticklabels():
#     tick.set_rotation(20)

print_title("Specificity of StruMs: auROC and auPRC")
t,p = stats.ttest_rel(x1,y1)
avg = np.average(x1-y1)
print """\
auROC::
    Avg Cross Family: {:>5.2f} +/- {:>5.2f}
    Avg Control:      {:>5.2f} +/- {:>5.2f}
    Avg improvement:  {:>5.2f} (t-stat: {:>5.2f}, p: {:>5.2e})
""".format(np.average(x1), np.std(x1), np.average(y1), np.std(y1), avg, t, p)
t,p = stats.ttest_rel(x2,y2)
avg = np.average(x2-y2)
print """\
auPRC:: 
    Avg Cross Family: {:>5.2f} +/- {:>5.2f}
    Avg Control:      {:>5.2f} +/- {:>5.2f}
    Avg improvement:  {:>5.2f} (t-stat: {:>5.2f}, p: {:>5.2e})
""".format(np.average(x2), np.std(x2), np.average(y2), np.std(y2), avg, t, p)

for ax in [ax_a, ax_b,]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
ax_b.set_yticklabels([])

plt.tight_layout()
fig3.savefig("figures/figure3.pdf")

######################################################################
# Figure 4: StruMs outperform, but are complementary to, PWMs (AUCs) #
######################################################################

fig4 = plt.figure(figsize=[twocol, 0.8*twocol])
plt.subplot(2,1,1)
label_plots(fig4, 2, 2)

pval_labels = ['P', 'D', 'S', 'C']

# 4A) auROCs of motifs using shuffled sequence as background
ax_a = plt.subplot(2,2,1)
plt.title("Shuffled background", weight='bold')
# plot_box(shuff_AUCs, pknonpk_labels)
plot_box(shuff_AUCs[:,[0,1,3,4]], pknonpk_labels[:2]+pknonpk_labels[3:])
plt.ylabel("auROC", weight='bold')

## STATISTICS
print_title("Difference in auROCs, shuffled background")
t_vals, p_vals = do_ttest(shuff_AUCs, pknonpk_labels2)
#plot_pvals(p_vals, [x[0] for x in pknonpk_labels], ax_a, cbar=False)
plot_pvals(p_vals[[0,1,3,4]][:,[0,1,3,4]], pval_labels, ax_a, cbar=False)

# 4B) auROCs of motifs using flanking sequence as background
ax_b = plt.subplot(2,2,2)
plt.title("Flanking background", weight='bold')
plot_box(flank_AUCs, pknonpk_labels)
# plt.ylabel("auROC", weight='bold')

## STATISTICS
print_title("Difference in auROCs, flanking seq as background")
t_vals, p_vals = do_ttest(flank_AUCs, pknonpk_labels2)
plot_pvals(p_vals, [x[0] for x in pknonpk_labels], ax_b)

# 4C) auPRCs of motifs using shuffled sequence as background
ax_c = plt.subplot(2,2,3)
plot_box(shuff_PRCs, pknonpk_labels)
plt.ylabel("auPRC", weight='bold')

## STATISTICS
print_title("Difference in auPRCs, shuffled background")
t_vals, p_vals = do_ttest(shuff_PRCs, pknonpk_labels2)
plot_pvals(p_vals, [x[0] for x in pknonpk_labels], ax_c, cbar=False)

# 4D) auPRCs of motifs using flanking sequence as background
ax_d = plt.subplot(2,2,4)
plot_box(flank_PRCs, pknonpk_labels)
# plt.ylabel("auPRC", weight='bold')

## STATISTICS
print_title("Difference in auPRCs, flanking seq as background")
t_vals, p_vals = do_ttest(flank_PRCs, pknonpk_labels2)
plot_pvals(p_vals, [x[0] for x in pknonpk_labels], ax_d)

plt.subplots_adjust(left=0.075, bottom=0.1, right=0.95, top=0.95,
                wspace=0.2, hspace=0.3)

for ax in [ax_a, ax_b, ax_c, ax_d,]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

plt.savefig("figures/figure4.pdf")

#######################################
# Figure 5: Complementarity of models #
#######################################

fig5 = plt.figure(figsize=[twocol, 2./3*twocol])
label_plots(fig5, 1, 3)

# 5A) Logit coefficients sorted by improvement over PWM AUC
pwm_auc = shuff_AUCs[:,0]
strum_auc = shuff_AUCs[:,3]

ax_a = plt.subplot(2,3,1)
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

plt.plot(log_pwm, log_strum, '.', c='grey', markersize=3)

plt.xlabel("Combined - PWM ($\Delta$ AUC)", weight='bold')
plt.ylabel("Combined - StruM ($\Delta$ AUC)", weight='bold')

plt.axhline(0, c='k', linewidth=1)
plt.axvline(0, c='k', linewidth=1)
ax_a.yaxis.set_label_coords(-0.175,0.5)

for side in ['top', 'right', 'left', 'bottom']:
        ax_a.spines[side].set_visible(False)


# 5B) Correlation of scores between sequence and shape models
R = []
for g in pknonpk_genes:
    if g in corr:
        R.append(corr[g])
left = 0.2
bottom = 0.16

x = R
y = log_pwm
ax_b = fig5.add_axes((0.33+0.33*left, 0.5+0.5*bottom, 0.33*(0.8-left), 0.5*(0.75-bottom)))
ax_pos = ax_b.get_position().bounds
ax_top = fig5.add_axes((0.33+0.33*left, 0.5+0.5*0.75, ax_pos[2], 0.5*0.2), sharex=ax_b)
ax_right = fig5.add_axes((0.33+0.33*0.8, 0.5+0.5*bottom, 0.33*0.2, ax_pos[3]), sharey=ax_b)
ax_b.plot(x, y, '.', markersize=3, c='grey')
ax_b.set_xlabel("Correlation between scores", weight='bold')
ax_b.set_ylabel("$\\Delta$ AUC (logit-PWM)", weight='bold')
ax_b.yaxis.set_label_coords(-0.175,0.5)
ax_top.hist(x, bins=30, facecolor='grey')
ax_right.hist(y, bins=30, facecolor='grey', orientation='horizontal')
ax_top.axis('off')
ax_right.axis('off')
z = np.polyfit(x, y, 1)
p = np.poly1d(z)
minx = min(x)
maxx = max(x)
ax_b.plot([minx,maxx],p([minx,maxx]), c='red', lw=1.5)

print_title("CORRELATION BETWEEN PWM AND STRUM SCORES")
print "Eqn of the line:   y = {:0.3f} x + {:0.3f}".format(z[0], z[1])
print "Correlation: {:0.3f} (p-value: {:0.2e})".format(*stats.pearsonr(x,y))
print "Avg. X (correlation):        {:>7.4f} +/- {:>5.2f}".format(np.average(x), np.std(x))
print "Avg. Y (Combined - PWM AUC): {:>7.4f} +/- {:>5.2f}".format(np.average(y), np.std(y))

# 5B) Avg distance from StruM to FIMO PWM
ax_c = plt.subplot(2,3,3)
plt.hist(avgdist_data, bins=np.linspace(0,60,31), normed=True)
ax_c.set_xlabel("Avg distance from top StruM\nmatch to sig PWM match (bp)", weight='bold')
ax_c.set_ylabel("Frequency", weight='bold', ha='right')
ax_c.yaxis.set_label_coords(-0.22, 0.75)

#5D) FIMO Positions - single peak
# Correlated
ax_d = plt.subplot(2,2,3)
fig5.text(0.01, 0.435, 'D', ha='left', va='top', weight='bold', fontsize=10)
cor_img = mpl.image.imread("output/" + "TAL1.ENCFF519DOC" # "eGFP-GTF2E2.ENCFF394WVZ" # "FOSL1.ENCFF961OPH"
                          + "_fimo_v_strum_matches2.png")
ax_d.set_title("TAL1", loc='left', fontsize=10)
ax_d.imshow(cor_img)
# ax_e.xaxis.set_visible(False)
ax_d.yaxis.set_visible(False)
xmin, xmax = ax_d.get_xlim()
xticks = range(int(xmin), int(xmax+1), int((xmax-xmin)/4))
xticks = np.linspace(xmin, xmax, 5)
ax_d.set_xticks(xticks)
ax_d.set_xticklabels(['-100', '0\nPWM Position (bp)', '+/-100', '0\nStruM Position (bp)', '+100'])

#5E) FIMO Positions - flanking peaks
# Flanking
ax_e = plt.subplot(2,2,4)
fig5.text(0.49, 0.435, 'E', ha='left', va='top', weight='bold', fontsize=10)
flank_img = mpl.image.imread("output/" + "RFX1.ENCFF934JXG" # "eGFP-ZNF512.ENCFF617CTX" # ZBTB33.ENCFF681IOP
                                + "_fimo_v_strum_matches2.png")
ax_e.set_title("RFX1", loc='left', fontsize=10)
ax_e.imshow(flank_img)
# ax_e.xaxis.set_visible(False)
ax_e.yaxis.set_visible(False)
xmin, xmax = ax_e.get_xlim()
xticks = range(int(xmin), int(xmax+1), int((xmax-xmin)/4))
xticks = np.linspace(xmin, xmax, 5)
ax_e.set_xticks(xticks)
ax_e.set_xticklabels(['-100', '0\nPWM Position (bp)', '+/-100', '0\nStruM Position (bp)', '+100'])


# # Anti Correlated
# ax_d = plt.subplot(3,1,4)
# anticor_img = mpl.image.imread("output/" + "MEF2A.ENCFF883WDT" # "RLF.ENCFF569QYK" # "ATF4.ENCFF491DNM" # ATF2.ENCFF525YRJ
#                           + "_fimo_v_strum_matches2.png")
# ax_d.set_title("MEF2A", loc='left')
# ax_d.imshow(anticor_img)
# ax_d.xaxis.set_visible(False)
# ax_d.yaxis.set_visible(False)

# # Not Correlated
# ax_e = plt.subplot(3,1,3)
# rand_img = mpl.image.imread("output/" + "SNIP1.ENCFF772GVZ" # "ATF4.ENCFF491DNM" # ATF2.ENCFF525YRJ
#                             + "_fimo_v_strum_matches2.png")
# ax_e.set_title("SNIP1", loc='left')
# ax_e.imshow(rand_img)
# # ax_e.xaxis.set_visible(False)
# ax_e.yaxis.set_visible(False)
# xmin, xmax = ax_e.get_xlim()
# xticks = range(int(xmin), int(xmax+1), int((xmax-xmin)/4))
# xticks = np.linspace(xmin, xmax, 5)
# ax_e.set_xticks(xticks)
# ax_e.set_xticklabels(['-100', '0\nPWM Position (bp)', '+/-100', '0\nStruM Position (bp)', '+100'])

# plt.subplots_adjust(left=0.1, bottom=0.075, right=0.9, top=0.975,
#                 wspace=0.2, hspace=0.3)

plt.subplots_adjust(left=0.075, bottom=0.05, right=0.95, top=0.95,
                wspace=0.3, hspace=0.3)

for ax in [ax_a, ax_c,]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

plt.savefig("figures/figure5.pdf", dpi=600)

#####################################
# Figure 6: Shape- vs. base-readout #
#####################################

# onecol = cm2inch(84/10.)
# twocol = cm2inch(178/10.)
# maxheight = cm2inch(230/10.)

oldheight = 1./3
newheight = 2./5

coeff_colors = ['#ff7f00', '#1f78b4']#['#ff7f00', '#33a02c', '#6a3d9a', '#1f78b4']
fig6 = plt.figure(figsize=[twocol, newheight*twocol])
# label_plots(fig6, 3, 1)

# 6A) Specificity of shape- vs base-readout TFs
ax_a = plt.subplot(1,3,1)

# goi_seq = ['STAT', 'GATA']
# goi_shp = ['TBP', 'LEF', 'RFX',]
goi_seq = [
    "STAT", "NFAT", "GATA", # Definitely in the dataset
    "BHTH", "homeodomains", # Classes of proteins
        "CTCF", #"CTCFL", #"DNMT1", "E4F1", "GATA2", "HINFP", "KLF13", # Zinc fingers
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

extra_colors = ["#fdbf6f", "#a6cee3"]#['gold', 'skyblue']
plt.plot(pwm_auc, strum_auc, '.', c=(169/256.,169/256.,169/256.), markersize=3,)
plt.plot(plt.xlim(), plt.ylim(), '--', c='grey')
plt.plot(pwm_auc[roi_seq], strum_auc[roi_seq], markers[1], 
         markerfacecolor=extra_colors[0], markeredgecolor='k', markeredgewidth=0.5, ms=4, label='Sequence')
plt.plot(pwm_auc[roi_shp], strum_auc[roi_shp], markers[2], 
         markerfacecolor=extra_colors[1], markeredgecolor='k', markeredgewidth=0.5, ms=4, label='Shape')


# for i,g,x,y,j in text:
#     plt.text(x, y, g, dict(weight='bold', size='x-small', ha=['left', 'right'][(j+i)%2], va=['top', 'bottom'][(j+i)%2]))

plt.legend(loc='upper left', fontsize='small', borderpad=0.2, labelspacing=0.2)
plt.xlabel("PWM auROC", weight='bold')
plt.ylabel("EM-StruM auROC", weight='bold')
ax_a.yaxis.set_label_coords(-0.175,0.5)
ax_a.set_aspect('equal', adjustable='box', anchor='C')

# 6B)  Logit coefficients sorted by improvement over PWM AUC


z = (newheight-oldheight)
left=0.075
# bottom=0.15
bottom= (z + 0.15*oldheight)/(newheight)
right=0.9
# top=0.9
top = (z + 0.9*oldheight)/(newheight)
wspace=0.3
hspace=0.3

fig6.text(0.01, top + (1-top)/2., 'A', ha='left', va='top', weight='bold', fontsize=10)
fig6.text(0.34, top + (1-top)/2., 'B', ha='left', va='top', weight='bold', fontsize=10)
ax_b = fig6.add_axes((0.37, bottom, right-0.37, top-bottom))
ax_below = fig6.add_axes((0.37, 0.0, right-0.37, bottom), sharex=ax_b)
ax_below.set_ylim((0,1))

ax_b.set_ylabel("logit coefficients", weight='bold',)
x = shuff_AUCs[:,4] - shuff_AUCs[:,0]
idx = np.argsort(x)

points = []
for i in [0,1]:
    p, = ax_b.plot(shuff_coeffs[:,i][idx], markers[[0,3][i]], c=coeff_colors[i], label=coeff_labels[i], markersize=3)
    points.append(p)

text = sorted(text, key=lambda x:(np.where(idx==x[0])[0]))


print_title("INFERRING READOUT MECHANISMS")
a = shuff_coeffs[:,0]
b = shuff_coeffs[:,1]
print "Number of experiments in which the logit model favors:"
print "(Any difference)"
print "- StruMs:  {}".format(sum(b > a))
print "- Neither: {}".format(sum(b == a))
print "- PWMs:    {}".format(sum(b < a))

print "(Two fold difference)"
print "- StruMs:  {}".format(sum((b / a) > 2)) ## TO DO: I NEED TO THINK ABOUT HOW THIS IS AFFECTED BY NEGATIVE VALUES
print "- Neither: {}".format(sum(((b / a) <= 2) & ((b / a) >= 0.5)))
print "- PWMs:    {}".format(sum((b / a) < 0.5))

minx, maxx = ax_b.get_xlim()
rangex = maxx-minx
c_width = 0.02*rangex
cluster_spacing = []

curr_clust = [ np.where(idx==text[0][0])[0][0], ]
curr_genes = [text[0][1],]
cmin = minx
i = 1
while i < len(text):
    ii = text[i][0]
    xval = np.where(idx==ii)[0][0]
    c_bottom = max(
                np.mean(curr_clust) - (c_width*2)*(len(curr_clust)/2.),
                cmin )
    if (xval) < (c_bottom + (c_width*2)*(len(curr_clust) + 0.5)):
        curr_clust.append(xval)
        curr_genes.append(text[i][1])
        i += 1
    else:
        for j,val in enumerate(curr_clust):
            cluster_spacing.append(c_bottom+j*(c_width*2))
        curr_clust = [ xval, ]
        curr_genes = [text[i][1]]
        cmin = c_bottom + (c_width*2)*(len(curr_clust) + 0.5)
        i += 1
c_bottom = max(
            np.mean(curr_clust) - (c_width*2)*(len(curr_clust)/2.),
            cmin )
for j,val in enumerate(curr_clust):
    cluster_spacing.append(c_bottom+j*(c_width*2))

for ii, (i,g,xval,yval,j) in enumerate(text):
    xval = np.where(idx==i)[0]
    y1 = shuff_coeffs[i,0]
    y2 = shuff_coeffs[i,1]
    yt = (y1+y2)/2.
    ys = [y1,y2]
    yt = ys[::-1][j]
    ax_b.plot([xval,xval], [y1,y2], '-', color=extra_colors[j], linewidth=1.5, 
                path_effects=[pe.Stroke(linewidth=2, foreground='k'), pe.Normal()])
    ax_b.plot([xval], [y1], 'd', ms=5, markerfacecolor=coeff_colors[0],
         markeredgecolor=extra_colors[j], path_effects=[pe.Stroke(linewidth=2, foreground='k'), pe.Normal()])
    ax_b.plot([xval], [y2], 'o', ms=5, markerfacecolor=coeff_colors[1],
         markeredgecolor=extra_colors[j], path_effects=[pe.Stroke(linewidth=2, foreground='k'), pe.Normal()])
    # ax_b.text(xval, yt, g, dict(weight='bold', size='small', ha='right',  va=['bottom', 'top', ][j], rotation=90))#ha=['left', 'right'][(j+i)%2], va=['top', 'bottom'][(j+i)%2]))
    c_pos = cluster_spacing[ii] + c_width
    ax_below.plot([xval, xval, c_pos, c_pos], [1.0, 0.95,0.85,0.81], 'k-', linewidth=0.5)
    ax_below.text(c_pos, 0.8, g, 
                    dict(weight='bold', size='small', ha='center',  
                        va='top', rotation=90,),
                    bbox=dict(facecolor=extra_colors[j], pad=1, edgecolor='None')
                    )
    # bbox=dict(facecolor='red', alpha=0.5)
ax_b2 = ax_b.twinx()
ax_b2.tick_params('y', colors='r')
ax_b2.set_ylabel("Combined - PWM ($\Delta$ AUC)", rotation=270, 
                 color='red', weight='bold', va='bottom',)
l, = plt.plot(x[idx], 'r-', label="Combined - PWM ($\Delta$ AUC)")
ax_b.set_xticks([])
ax_b.set_xlabel("Ranked by $\Delta$ AUC", weight='bold')
ax_b2.legend(points, [coeff_labels[i] for i in [0,1]], bbox_to_anchor=[0.5, 1.0], loc='center', ncol=2, framealpha=1, facecolor='white')

ax_b.yaxis.set_label_coords(-0.045,0.5)
ax_b2.yaxis.set_label_coords(1.1,0.5)
ax_b.xaxis.set_label_coords(0.5,0.08)
# ax_a.xaxis.set_label_coords(0.5,0.08)

ax_below.spines["top"].set_visible(False)
ax_below.spines["right"].set_visible(False)
ax_below.spines["left"].set_visible(False)
ax_below.spines["bottom"].set_visible(False)
ax_below.set_yticks([])

#### Adjust Subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top,
                wspace=wspace, hspace=hspace)

plt.savefig('figures/figure6.pdf')

######################################################


dist_data = dict(zip(avgdist_tfs, avgdist_data))
x = [dist_data[g] for g in pknonpk_genes if g in dist_data]
y = [shuff_AUCs[i,3] for i,g in enumerate(pknonpk_genes) if g in dist_data]



sfig = plt.figure(figsize=(onecol, onecol))
plt.plot(x, y, '.', color=(169/256.,169/256.,169/256.))

z = np.polyfit(x, y, 1)
p = np.poly1d(z)
minx = min(x)
maxx = max(x)

plt.plot([minx,maxx],p([minx,maxx]), c='red', lw=1.5)
plt.xlabel("Avg distance between best StruM score\nand nearest FIMO position")
plt.ylabel("EM-StruM auROC (shuff)")
plt.tight_layout()

r,p = stats.pearsonr(x,y)

print_title("RELATION BETWEEN DISTANCE FROM FIMO-PWM SITE AND EM-STRUM PERFORMANCE")
print "Eqn of the line:   y = {:0.3f} x + {:0.3f}".format(z[0], z[1])
print "Correlation: {:0.3f} (p-value: {:0.2e})".format(r, p)
print "Avg. X (Avg Dist):               {:>7.4f} +/- {:>5.2f}".format(np.average(x), np.std(x))
print "Avg. Y (EM-StruM auROC (shuff)): {:>7.4f} +/- {:>5.2f}".format(np.average(y), np.std(y))

plt.savefig('Supplemental/fimo_perf.pdf')
