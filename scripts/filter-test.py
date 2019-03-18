import numpy as np
import matplotlib.pyplot as plt
import string
import sys

from scipy.special import ndtr as norm_p
from scipy.interpolate import UnivariateSpline
from statsmodels.nonparametric.kde import KDEUnivariate
from sklearn.metrics import precision_recall_curve, roc_curve, auc

from strum import strum

##################################################
# DEFINE HELPER FUNCTIONS
##################################################
def fisher_score(x0, x1):
    mu0, mu1 = np.average(x0), np.average(x1)
    sig0, sig1 = np.std(x0), np.std(x1)
    return (mu0-mu1)**2/(sig0**2+sig1**2)

def filt_eval(struc_kmer, min_i, idx):
    filt = idx[min_i:]
    return np.sum(np.log10(norm_p(
        -np.absolute(struc_kmer[filt]-motif.strum[0][filt])/motif.strum[1][filt]
    ) + 10**-300))

def scale_0to1(x):
    x -= np.min(x)
    x /= np.max(x)
    return x

def cm2inch(value):
    return value/2.54

def motif_eval(x, mu, sig):
    return np.sum(np.log10(norm_p(-np.absolute(x-mu)/sig)+ 10**-300)) 

onecol = cm2inch(84/10.)
twocol = cm2inch(178/10.)
maxheight = cm2inch(230/10.)

def print_title(string):
    l = len(string) + 4
    print "\n" + "#"*l
    print "# {} #".format(string.upper())
    print "#"*l

##################################################
# Train the motif and score sequences
##################################################
k = 11
motif = strum.StruM(mode='proteingroove')

tf_file = sys.argv[1]
print tf_file

with open(tf_file) as f:
    header, sequences = motif.read_FASTA(f)
top_seq = [x.upper() for x in sequences[:500]]
next_seq = [x.upper() for x in sequences[500:]]
motif.train_EM(top_seq, fasta=False, lim=0.1, 
        k=(k-1), max_iter=250, n_init=5, random_seed=5407)
n = len(motif.strum[0])

sequences = []
for seq in top_seq:
    rseq = strum.rev_comp(seq)
    fscore = motif.score_seq(seq)
    rscore = motif.score_seq(rseq)[::-1]
    fi = np.argmax(fscore)
    ri = np.argmax(rscore)
    if fscore[fi] > rscore[ri]:
        sequences.append(seq[fi:fi+k])
    else:
        sequences.append(rseq[ri:ri+k])

neg_seqs = []
for seq in sequences:
    tmp_seq = list(seq)
    np.random.shuffle(tmp_seq)
    neg_seqs.append(''.join(tmp_seq))

motif.print_PWM()


##################################################
# COMPUTE THE FISHER SCORE FOR EACH POSITION-
# SPECIFIC FEATURE
##################################################

pos_shape = np.array([motif.translate(x) for x in sequences])
neg_shape = np.array([motif.translate(x) for x in neg_seqs])

V = []
for i in range(pos_shape.shape[1]):
    V.append(fisher_score(pos_shape[:,i], neg_shape[:,i]))
V = np.array(V)

# Sort 
idx = np.argsort(V)
sortV = V[idx]
logSortV = np.log(sortV)

# Fit a spline to automatically choose a threshold

xvals = range(len(logSortV))
spl = UnivariateSpline(xvals, logSortV, s=n/10.)
ys = spl(xvals)

d_spl = spl.derivative(1)
d_ys = d_spl(xvals)

min_i = np.argmin(d_ys)

plt.figure(figsize=(twocol,twocol/2.))
plt.subplot(121)
plt.title("Fit")
plt.plot(logSortV)
plt.plot(ys)
plt.axvline(min_i, ls='--', c='gray')
plt.xlabel("Position-specific Feature Rank")
plt.ylabel("Log Fisher Score")
plt.subplot(122)
plt.title("Residuals")
plt.plot(xvals, ys-logSortV, '.')
plt.axhline(0, ls='--', c='gray');
plt.xlabel("Position-specific Feature Rank")
plt.ylabel("Residual Log Fisher Score")
plt.tight_layout()
plt.savefig('Supplemental/filt1.pdf')
# print(spl(min_i))

##################################################
# CHOOSE A VARIANCE THRESHOLD TO FILTER ON
##################################################

idx2 = np.argsort(motif.strum[1])[::-1]
variance = motif.strum[1][idx2]

spl2 = UnivariateSpline(xvals, variance, s=n/10.)
ys2 = spl2(xvals)

d_spl2 = spl2.derivative(1)
d_ys2 = d_spl2(xvals)

min_i2 = np.argmax(d_ys2)

plt.figure(figsize=(twocol,twocol/2.))
plt.subplot(121)
plt.title("Fit")
plt.plot(variance)
plt.plot(ys2)
plt.axvline(min_i2, ls='--', c='gray')
plt.xlabel("Position-specific Feature Rank")
plt.ylabel("Standard Deviation")
plt.subplot(122)
plt.title("Residuals")
plt.plot(xvals, ys2-variance, '.')
plt.axhline(0, ls='--', c='gray');
plt.xlabel("Position-specific Feature Rank")
plt.ylabel("Residual Standard Deviation")
plt.tight_layout()
plt.savefig('Supplemental/filt2.pdf')
# print(spl2(min_i2))

##################################################
# COMPARE FILTERING SCHEMES
##################################################

Y = np.array([1]*len(pos_shape) + [0]*len(neg_shape))
full = np.array([motif_eval(x, motif.strum[0], motif.strum[1]) for x in pos_shape] + [motif_eval(x, motif.strum[0], motif.strum[1]) for x in neg_shape])
filtered = np.array([filt_eval(x, min_i, idx) for x in pos_shape] + [filt_eval(x, min_i, idx) for x in neg_shape])
filtered2 = np.array([filt_eval(x, min_i2, idx2) for x in pos_shape] + [filt_eval(x, min_i2, idx2) for x in neg_shape])

### SEPARABILITY ###

print_title("SEPARABILITY")

print "| {:<31s} | {:6} |".format("Version", "Fisher")
print "|:{:<31s}-|:{:6}:|".format("-"*31, "-"*5,)
print "| {:<31s} | {:6.3f} |".format("ProteinGroove", fisher_score(full[Y==0], full[Y==1]))
print "| {:<31s} | {:6.3f} |".format("ProteinGroove - Fisher Filter", fisher_score(filtered[Y==0], filtered[Y==1]))
print "| {:<31s} | {:6.3f} |".format("ProteinGroove - Variance Filter", fisher_score(filtered2[Y==0], filtered2[Y==1]))


### AUC ###


full_0to1 = scale_0to1(full)
filtered_0to1 = scale_0to1(filtered)
filtered2_0to1 = scale_0to1(filtered2)

full_precision, full_recall, thresholds = precision_recall_curve(Y, full_0to1)
full_fpr, full_tpr, thresholds = roc_curve(Y, full_0to1)
filtered_precision, filtered_recall, thresholds = precision_recall_curve(Y, filtered_0to1)
filtered_fpr, filtered_tpr, thresholds = roc_curve(Y, filtered_0to1)
filtered2_precision, filtered2_recall, thresholds = precision_recall_curve(Y, filtered2_0to1)
filtered2_fpr, filtered2_tpr, thresholds = roc_curve(Y, filtered2_0to1)

fig = plt.figure(figsize=(twocol,twocol/2.))
ax1 = plt.subplot(121)
ax1.plot(full_fpr, full_tpr, label='PG')
ax1.plot(filtered_fpr, filtered_tpr, label='PG-Fish')
ax1.plot(filtered2_fpr, filtered2_tpr, label='PG-Var')
ax1.set_xlabel("FPR")
ax1.set_ylabel("TPR")
ax2 = plt.subplot(122)
ax2.plot(full_recall, full_precision, label='PG')
ax2.plot(filtered_recall, filtered_precision, label='PG-Fish')
ax2.plot(filtered2_recall, filtered2_precision, label='PG-Var')
ax2.set_xlabel("Recall")
ax2.set_ylabel("Precision")
ax2.legend();
plt.tight_layout()

print_title("PREDICTIVE PERFORMANCE")

print "| {:<31s} | {:6} | {:5} |".format("Version", "auROCC", "auPRC")
print "|:{:<31s}-|:{:6}:|:{:5}:|".format("-"*31, "-"*6, "-"*5)
print "| {:<31s} | {:6.3f} | {:5.3f} |".format("ProteinGroove", auc(full_fpr, full_tpr), auc(full_recall, full_precision))
print "| {:<31s} | {:6.3f} | {:5.3f} |".format("ProteinGroove - Fisher Filter", auc(filtered_fpr, filtered_tpr), auc(filtered_recall, filtered_precision))
print "| {:<31s} | {:6.3f} | {:5.3f} |".format("ProteinGroove - Variance Filter", auc(filtered2_fpr, filtered2_tpr), auc(filtered2_recall, filtered2_precision))


##################################################
# COMPARE MODES FULL v PROTEINGROOVE
##################################################



motif = strum.StruM(mode='full')
motif.train(sequences, fasta=False)
n = len(motif.strum[0])

idx2 = np.argsort(motif.strum[1])[::-1]
variance = motif.strum[1][idx2]

xvals = range(len(variance))
spl2 = UnivariateSpline(xvals, variance, s=n/10.)
ys2 = spl2(xvals)

d_spl2 = spl2.derivative(1)
d_ys2 = d_spl2(xvals)

min_i2 = np.argmax(d_ys2)

pos_shape = np.array([motif.translate(x) for x in sequences])
neg_shape = np.array([motif.translate(x) for x in neg_seqs])

full = np.array([motif_eval(x, motif.strum[0], motif.strum[1]) for x in pos_shape] + [motif_eval(x, motif.strum[0], motif.strum[1]) for x in neg_shape])
filtered3 = np.array([filt_eval(x, min_i2, idx2) for x in pos_shape] + [filt_eval(x, min_i2, idx2) for x in neg_shape])
full_0to1 = scale_0to1(full)
filtered3_0to1 = scale_0to1(filtered3)
full_precision, full_recall, thresholds = precision_recall_curve(Y, full_0to1)
full_fpr, full_tpr, thresholds = roc_curve(Y, full_0to1)
filtered3_precision, filtered3_recall, thresholds = precision_recall_curve(Y, filtered3_0to1)
filtered3_fpr, filtered3_tpr, thresholds = roc_curve(Y, filtered3_0to1)

print "| {:<31s} | {:6.3f} | {:5.3f} |".format("Full", auc(full_fpr, full_tpr), auc(full_recall, full_precision))
ax1.plot(full_fpr, full_tpr, '--', label='Full')
ax2.plot(full_recall, full_precision, '--', label='Full')

print "| {:<31s} | {:6.3f} | {:5.3f} |".format("Full - Variance Filter", auc(filtered3_fpr, filtered3_tpr), auc(filtered3_recall, filtered3_precision))
ax1.plot(filtered3_fpr, filtered3_tpr, '--', label='Full-Var')
ax2.plot(filtered3_recall, filtered3_precision, '--', label='Full-Var')

ax2.legend();
fig.savefig('Supplemental/filt3.pdf')