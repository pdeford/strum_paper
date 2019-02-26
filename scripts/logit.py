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
from sklearn.metrics import roc_curve, auc, precision_recall_curve
# from sklearn.cross_validation import StratifiedKFold # Version 0.17.1
from sklearn.model_selection import StratifiedKFold    # Version 0.20.1

try:
    import cPickle as pickle
except:
    import pickle

from strum import strum


# Set up indexes for defining the weight matrices
nucs = "ACGT"
nuc_index = dict(zip(nucs, range(4)))
dimers = [a+b for a in nucs for b in nucs]
di_index = dict(zip(dimers, range(16)))
N_seq = 500

def score_all(basename, n_process, random_seed, models, sequences):
    global k, pwm, dwm, ml_strum, em_strum

    random.seed(random_seed)

    pwm, dwm, ml_strum, em_strum = models
    k = pwm.shape[1]
    N_seq = 500
    print >> sys.stderr, "Load positive and negative sequences"
    # Get positive examples from file
    peaks = sequences[N_seq:2*N_seq]

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

    data2 = []
    for i, seq_set in enumerate([decoys, peaks]):
        data_ = []
        for seq in seq_set:
            ml_strum_scores = [
                ml_strum.score_seq_filt(seq),
                ml_strum.score_seq_filt(rev_comp(seq))
                ]
            em_strum_scores = [
                em_strum.score_seq_filt(seq),
                em_strum.score_seq_filt(rev_comp(seq))
                ]
            data_.append((np.max(np.hstack(ml_strum_scores)), np.max(np.hstack(em_strum_scores))))
        data2 += data_



    data = np.vstack(data)  
    data2 = np.vstack(data2)
    data = np.hstack([data, data2])
    y = np.array(y, dtype=int)

    # Filter out bad sequences
    bad_rows = np.where(np.isnan(data) | np.isinf(data))[0]
    mask = np.ones(data.shape[0], dtype=bool)
    mask[bad_rows] = False
    data = data[mask]
    y = y[mask]
    print >> sys.stderr, "Skipping {} sequences due to `nan` or `inf`".format(np.sum(~mask))

    print >> sys.stderr, "Get ROC curves"
    # Evaluate the performance of each motif
    results = []
    for i, motif in enumerate([pwm, dwm, ml_strum, em_strum]):
        fpr, tpr, _ = roc_curve(y, data[:,i])
        precision, recall, _ = precision_recall_curve(y, data[:,i])
        results.append((fpr, tpr, auc(fpr,tpr), recall, precision, auc(recall, precision)))

    ############## REPEAT FOR OTHER TYPES OF DECOYS #############
    print >> sys.stderr, "Load other decoys"
    decoys2 = fasta_reader(open("data/{}.flank.fa".format(basename)))[N_seq:]

    print >> sys.stderr, "Score decoys"
    pool = Pool(n_process)
    data2 = pool.map(eval_seq, decoys2)
    pool.close()
    pool.join()
    data2 = np.vstack(data2)

    data3 = []
    for seq in decoys2:
        ml_strum_scores = [
            ml_strum.score_seq_filt(seq),
            ml_strum.score_seq_filt(rev_comp(seq))
            ]
        em_strum_scores = [
            em_strum.score_seq_filt(seq),
            em_strum.score_seq_filt(rev_comp(seq))
            ]
        data3.append((np.max(np.hstack(ml_strum_scores)), np.max(np.hstack(em_strum_scores))))



    data = np.vstack(data)  
    data2 = np.vstack(data2)
    data3 = np.vstack(data3)
    data2 = np.hstack([data2, data3])
    data = np.hstack([data, data2])

    print >> sys.stderr, "Decoys ROCs"
    results2 = []
    working_x = np.vstack([data[y == 1], data2])
    working_y = np.zeros([len(working_x),], dtype=int)
    working_y[:sum(y)] = 1
    for i, motif in enumerate([pwm, dwm, ml_strum, em_strum]):
        fpr, tpr, _ = roc_curve(working_y, working_x[:,i])
        precision, recall, _ = precision_recall_curve(working_y, working_x[:,i])
        results2.append((fpr, tpr, auc(fpr,tpr), recall, precision, auc(recall, precision)))

    ##############################################################

    print >> sys.stderr, "Logit"
    # Train a combined model using logistic regression
    ## Scale the input so each set of scores ranges in [0,1]
    scaler = (np.min(data, axis=0), np.max(data, axis=0))
    data -= scaler[0]
    data /= (scaler[1] - scaler[0])

    ## Train the model
    clf = logit()
    # cv = StratifiedKFold(y,n_folds=10)  # Version 0.17.1
    cv = StratifiedKFold(n_splits=10)     # Version 0.20.1
    tprs = []
    precisions = []
    aucs = []
    prcaucs = []
    mean_fpr = np.linspace(0, 1, 100)
    # for train, test in cv:              # Version 0.17.1
    for train, test in cv.split(data, y): # Version 0.20.1
        probas_ = clf.fit(data[train], y[train]).predict_proba(data[test])
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
        precision, recall, thresholds = precision_recall_curve(y[test], probas_[:, 1])
        tprs.append(interp(mean_fpr, fpr, tpr))
        precisions.append(interp(mean_fpr, recall, precision))
        tprs[-1][0] = 0.0
        precisions[-1][-1] = 0.0
        roc_auc = auc(fpr, tpr)
        prcauc = auc(recall, precision)
        aucs.append(roc_auc)
        prcaucs.append(prcauc)
    mean_tpr = np.mean(tprs, axis=0)
    mean_precision = np.mean(precisions, axis=0)
    mean_tpr[-1] = 1.0
    mean_precision[0] = 1.0
    results.append((mean_fpr, mean_tpr, np.mean(aucs), mean_fpr, mean_precision, np.mean(prcaucs))) # Combined model


    ############## REPEAT FOR OTHER TYPES OF DECOYS #############
    # Train a combined model using logistic regression
    ## Scale the input so each set of scores ranges in [0,1]
    print >> sys.stderr, "...and Again"
    scaler2 = (np.min(working_x, axis=0), np.max(working_x, axis=0))
    working_x -= scaler2[0]
    working_x /= (scaler2[1] - scaler2[0])

    ## Train the model
    clf2 = logit()
    # cv2 = StratifiedKFold(working_y, n_folds=10)      # Version 0.17.1
    cv2 = StratifiedKFold(n_splits=10)                  # Version 0.20.1

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    # for train, test in cv2:                           # Version 0.17.1
    for train, test in cv2.split(working_x, working_y): # Version 0.20.1
        probas_ = clf2.fit(working_x[train], working_y[train]).predict_proba(working_x[test])
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(working_y[test], probas_[:, 1])
        precision, recall, thresholds = precision_recall_curve(working_y[test], probas_[:, 1])
        tprs.append(interp(mean_fpr, fpr, tpr))
        precisions.append(interp(mean_fpr, recall, precision))
        tprs[-1][0] = 0.0
        precisions[-1][-1] = 0.0
        roc_auc = auc(fpr, tpr)
        prcauc = auc(recall, precision)
        aucs.append(roc_auc)
        prcaucs.append(prcauc)
    mean_tpr = np.mean(tprs, axis=0)
    mean_precision = np.mean(precisions, axis=0)

    mean_tpr[-1] = 1.0
    mean_precision[0] = 1.0
    results2.append((mean_fpr, mean_tpr, np.mean(aucs), mean_fpr, mean_precision, np.mean(prcaucs))) # Combined model

    ##############################################################

    print >> sys.stderr, "Print Results"
    # FORMAT:
    ##FACTOR.ACCESSION  pwm-auROC   dwm-auROC   ml_strum-auROC  em_strum-auROC  logit-auROC pwm-auPRC   dwm-auPRC   ml_strum-auPRC  em_strum-auPRC  logit-auPRC
    # Output results and plot the results
    print basename + "\t" + "\t".join(["%0.6f" % x[2] for x in results]) + "\t" + "\t".join(["%0.6f" % x[5] for x in results])
    ############## REPEAT FOR OTHER TYPES OF DECOYS #############
    print basename + "\t" + "\t".join(["%0.6f" % x[-1] for x in results2]) + "\t" + "\t".join(["%0.6f" % x[5] for x in results2])
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

    plt.figure()
    for i, name in enumerate(labels):
        row = results[i]
        plt.plot(row[3], row[4], c=colors[i], label=name)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    #plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.savefig("output/{}_PRC.pdf".format(basename))
    plt.close()

    plt.figure()
    for i, name in enumerate(labels):
        row = results[i]
        plt.plot(row[3], row[4], c=colors[i], label=name)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    #plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.savefig("output/{}_PRC2.pdf".format(basename))
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
        for j in range(len(seq) - k + 1):
            kmer = seq[j:j+k]
            if 'N' in kmer: continue
            rkmer = rev_comp(kmer)
            pwm_scores.append(score_pwm(pwm, kmer))
            pwm_scores.append(score_pwm(pwm, rkmer))
            dwm_scores.append(score_dwm(dwm, kmer))
            dwm_scores.append(score_dwm(dwm, rkmer))
        return (np.max(pwm_scores), np.max(dwm_scores))

if __name__ == '__main__':
    basename = sys.argv[1]
    n_process = int(sys.argv[2])
    random_seed = int(sys.argv[3])
    models = pickle.load(open(sys.argv[4], 'rb'))
    
    score_all(basename, n_process, random_seed, models)
