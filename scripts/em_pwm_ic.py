#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')

import os
import numpy as np
from scipy.stats import pearsonr

import sys

try:
    import cPickle as pickle
except:
    import pickle

import matplotlib.pyplot as plt
import matplotlib.patches as patches

line_text = """\
y = {:0.2f} x + {:0.2f}
R: {:0.2f}"""

def main():
    aucs = {}
    keepers = [0,3,4]
    auc_file = open(sys.argv[1])
    i = 0
    for line in auc_file:
        fields = line.split()
        tf = fields[0]
        auc = [float(x) for x in fields[1:]]
        if i % 2 == 0:
            aucs[tf] = [auc[x] for x in keepers]
        # else:
        #     aucs[tf] += [auc[c] for x in keepers]
    motifs = []
    for motif_file in sys.argv[2:]:
        with open(motif_file, "rb") as f:
            pwm, dwm, ml_strum, em_strum, (logit, scaler), (logit2, scaler2) = pickle.load(f)
        mname = motif_file.split("/")[-1].rstrip('.p')

        rpwm = rev_comp_pwm(pwm)
        i, alignment = pwm_similarity(pwm, em_strum.PWM)
        i2, alignment2 = pwm_similarity(rpwm, em_strum.PWM)
        if alignment2 > alignment:
            pwm, i, alignment = rpwm, i2, alignment2
        p_ic = pwm_ic(pwm)
        s_ic = pwm_ic(em_strum.PWM)
        k = pwm.shape[1]

        if mname in aucs:
            motifs.append(
                [mname, pwm, em_strum.PWM, k, i, alignment, p_ic, np.average(p_ic), np.sum(p_ic), s_ic, np.average(s_ic), np.sum(s_ic)] + aucs[mname]
            )

    motifs.sort(key=lambda x:x[5:8])
    motifs = motifs[::-1]

    f1 = plt.figure(figsize=[10,5])
    x = [tup[7] for tup in motifs]
    y = [tup[8] for tup in motifs]
    z = [tup[3] for tup in motifs]
    ax1 = f1.add_subplot(121)
    scat = ax1.scatter(x,y,c=z, cmap='viridis')
    ax1.set_xlabel("Average IC")
    ax1.set_ylabel("Total IC")
    ax1.set_title('PWM')
    x = [tup[10] for tup in motifs]
    y = [tup[11] for tup in motifs]
    z = [tup[3] for tup in motifs]
    ax2 = f1.add_subplot(122)
    scat = ax2.scatter(x,y,c=z, cmap='viridis')
    ax2.set_xlabel("Average IC")
    ax2.set_ylabel("Total IC")
    ax2.set_title('StruM.PWM')
    cbar = plt.colorbar(scat, ax=ax2)
    cbar.set_label('Motif width', rotation=270, va='bottom')
    f1.savefig('output/pwm_qual.png')

    f2 = plt.figure()
    ax3 = f2.add_subplot(111)
    x = [tup[5] for tup in motifs]
    ax3.hist(x, bins=30)
    f2.savefig('output/motif_similarity_dist.png')

    nucs = [plot_A, plot_C, plot_G, plot_T]
    try:
        os.mkdir(os.path.join('output','motif_alignments'))
    except:
        pass
    for mname, pwm, s_pwm, k, a_i, alignment, p_ic, p_ic_avg, p_ic_total, s_ic, s_ic_avg, s_ic_total, pauroc, sauroc, lauroc in motifs:
        f = plt.figure(figsize=[k/2.,(2*4)/2.])
        ax_p = f.add_subplot(211)
        for x in range(k):
            y = 0
            idx = np.argsort(pwm[:,x],)
            for i in idx:
                height = pwm[i,x]*p_ic[x]
                nucs[i](ax_p, x, y, height)
                y += height
        ax_p.set_ylabel('PWM IC (bits)')
        ax_p.set_xticks(np.arange(k)+0.5)
        ax_p.set_xticklabels(range(k))
        ax_s = f.add_subplot(212)
        for x in range(k):
            y = 0
            idx = np.argsort(s_pwm[:,x],)
            for i in idx:
                height = s_pwm[i,x]*s_ic[x]
                nucs[i](ax_s, x+a_i, y, height)
                y += height
        ax_s.set_ylabel('StruM.PWM IC (bits)')
        ax_s.set_xticks(np.arange(k)+0.5+a_i)
        ax_s.set_xticklabels(range(k))
        minx = min(0, a_i)
        maxx = max(k, a_i+k)
        ax_p.set_xlim((minx, maxx))
        ax_p.set_ylim((0,2))
        ax_s.set_xlim((minx, maxx))
        ax_s.set_ylim((0,2))
        f.suptitle("{}\nScore: {:0.3f}".format(mname, alignment))
        # f.suptitle(mname)

        f.savefig(os.path.join('output', 'motif_alignments', '{}.png'.format(mname)))
        plt.close(f)

    n = len(motifs)
    a = int(np.sqrt(n)//1)
    b = int(n//a + (n%a > 0))
    if a*b < n:
        b += 1

    diff = a*b - n
    datums = {}
    #0      1     2    3   4     5         6      7          8         9       10        11          12     13      14
    #mname, pwm, s_pwm, k, a_i, alignment, p_ic, p_ic_avg, p_ic_total, s_ic, s_ic_avg, s_ic_total, pauroc, sauroc, lauroc
    labels = ['Alignment Score', 'PWM Avg IC', 'PWM Total IC', 'Motif Width', 'StruM.PWM Avg IC', 'StruM.PWM Total IC']
    indices = [5, 7, 8, 3, 10, 11]
    labels = ['PWM Avg IC', 'StruM.PWM Avg IC', 'Alignment Score', 
              'PWM Total IC', 'StruM.PWM Total IC', 'Motif Width',
              'PWM auROC', 'StruM auROC', 'logit auROC',
              #'PWM auPRC', 'StruM auPRC', 'logit auPRC',
              ]
    indices = [7,  10, 5,
               8,  11, 3,
               12, 13, 14,
               #15, 16, 17,
               ]
    for i,l in zip(indices, labels):
        v = [m[i] for m in motifs]
        for _ in range(diff): v.append(np.nan)
        datums[l] = v
    f3 = plt.figure(figsize=(4*a,2*b))
    for i in range(n):
        mname, pwm, s_pwm, k, a_i, alignment, p_ic, p_ic_avg, p_ic_total, s_ic, s_ic_avg, s_ic_total, pauroc, sauroc, lauroc = motifs[i]
        ax = f3.add_subplot(b, a, i+1)
        ## Generate plot of aligned motifs
        for x in range(k):
            y = 0
            idx = np.argsort(pwm[:,x],)
            for i in idx:
                height = pwm[i,x]*p_ic[x]
                nucs[i](ax, x, y, height)
                y += height
        for x in range(k):
            y = -2
            idx = np.argsort(s_pwm[:,x],)
            for i in idx:
                height = s_pwm[i,x]*s_ic[x]
                nucs[i](ax, x+a_i, y, height)
                y += height
        ax.set_ylim((-2,2))
        minx = min(0, a_i)
        maxx = max(k, a_i+k)
        ax.set_xlim((minx, maxx))
        ax.set_xlim()
        ax.set_xticks(())
        ax.set_yticks(())
        ax.set_title(mname)
    plt.tight_layout()
    f3.savefig('output/sorted_motifs.png')

    
    f4 = plt.figure(figsize=(12,6))
    ncol = 4
    nrow = 3
    for i, l in enumerate(labels):
        idx = int(i + i//(ncol-1))
        ax = plt.subplot(nrow,ncol,idx+1)
        plt.title(l)
        plt.pcolor(np.reshape(datums[l], (b,a)))
        ax.invert_yaxis()
        plt.colorbar()
    ax = plt.subplot(1,ncol,ncol)
    big_arr = np.vstack([datums[l][:n] for l in labels]).T
    big_arr -= np.min(big_arr, axis=0)
    big_arr /= np.max(big_arr, axis=0)
    plt.pcolor(big_arr)
    ax.invert_yaxis()
    plt.tight_layout()
    f4.savefig('output/sorted_motif_heatmaps.png')

    m = len(labels)
    f5 = plt.figure(figsize=(3*m,3*m))
    f5, axes = plt.subplots(nrows=m, ncols=m, figsize=(3*m,3*m))
    for i in range(m):
        y = datums[labels[i]][:n]
        for j in range(m):
            x = datums[labels[j]][:n]
            ax = axes[i,j]
            if i == j:
                ax.hist(y, bins=11, facecolor='#4E73AE', edgecolor='white')
                if j > 0:
                    ax.get_shared_x_axes().join(ax, axes[0,j])
            else:
                ax.plot(x,y, '.', c='#4E73AE')
                z = np.polyfit(x, y, 1)
                r = np.corrcoef(x,y)[0,1]
                ax.text(1, 1, line_text.format(z[0], z[1], r), 
                              ha='right', va='top', transform=ax.transAxes,
                              bbox=dict(facecolor='white', edgecolor='white', alpha=0.5))
                if i > 0:
                    ax.get_shared_x_axes().join(ax, axes[0,j])
                if j > 0:
                    if i == 0:
                        if j > 1:
                            ax.get_shared_y_axes().join(ax, axes[i,1])
                    else:
                        ax.get_shared_y_axes().join(ax, axes[i,0])
            ax.autoscale()
            for side in ['top', 'right',]:
                ax.spines[side].set_visible(False)
            if i < m-1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel(labels[j], weight='bold')
            if j > 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel(labels[i], weight='bold')
    plt.tight_layout()
    f5.savefig('output/scatter_matrix.png')

def pwm_similarity(m1, m2):
    k1 = m1.shape[1]
    k2 = m2.shape[1]

    assert (k1 == k2)

    compiled_d = [-1]
    for start in range(-k2 + 1, k1):
        sub_d = 0.
        for i in range(k2):
            if ((i + start) >= 0) and ((start + i) < k1):
                r, p = pearsonr(m2[:,i], m1[:, start+i])
                if np.isnan(r):
                    break
                sub_d += r
        else:
            compiled_d.append(sub_d)

    i = np.argmax(compiled_d)
    # Position relative to start of m1 corresponding to first position of m2
    # Also the score of that alignment
    return -k2 + i, compiled_d[i]

def rev_comp_pwm(pwm):
    return pwm[::-1, ::-1]

def pwm_ic(pwm):
    # Add pseudo-counts, just in case
    pwm += 0.00001
    pwm /= np.sum(pwm, axis=0)

    # Find information content at each position
    ic = np.log2(pwm)* pwm
    ic = np.sum(ic, axis=0) + 2.0

    return ic

def plot_A(ax, x, y, height, width=1):
    xs = [0.0, 0.5, 1.0, 0.85, 0.65, 0.35, 0.15, 0.0]
    ys = [0.0, 1.0, 0.0, 0.0, 0.4, 0.4, 0.0, 0.0]

    xs = [i*width+x for i in xs]
    ys = [i*height+y for i in ys]
    ax.add_patch(patches.Polygon(xy=list(zip(xs,ys)), fill=True, facecolor='green'))

    xs = [0.4, 0.5, 0.6, 0.4]
    ys = [0.55, 0.8, 0.55, 0.55]
    xs = [i*width+x for i in xs]
    ys = [i*height+y for i in ys]
    ax.add_patch(patches.Polygon(xy=list(zip(xs,ys)), fill=True, facecolor='white'))

def plot_C(ax, x, y, height, width=1):
    cx = x + width/2.
    cy = y + height/2.

    ax.add_patch(patches.Ellipse((cx,cy), width, height, facecolor='blue'))
    ax.add_patch(patches.Ellipse((cx,cy), width*0.7, height*0.7, facecolor='white'))
    ax.add_patch(patches.Rectangle((cx,cy-0.15*height), width*0.5, 0.3*height, facecolor='white'))

def plot_G(ax, x, y, height, width=1):
    cx = x + width/2.
    cy = y + height/2.

    ax.add_patch(patches.Ellipse((cx,cy), width, height, facecolor='gold'))
    ax.add_patch(patches.Ellipse((cx,cy), width*0.7, height*0.7, facecolor='white'))
    ax.add_patch(patches.Rectangle((cx,cy-0.15*height), width*0.5, 0.3*height, facecolor='white'))
    ax.add_patch(patches.Rectangle((cx,cy-0.15*height), width*0.5, 0.15*height, facecolor='gold'))

def plot_T(ax, x, y, height, width=1):
    xs = [0, 1, 1, 0.6, 0.6, 0.4, 0.4, 0, 0 ]
    ys = [1, 1, 0.8, 0.8, 0, 0, 0.8, 0.8, 1]

    xs = [i*width+x for i in xs]
    ys = [i*height+y for i in ys]
    ax.add_patch(patches.Polygon(xy=list(zip(xs,ys)), fill=True, facecolor='red'))

if __name__ == '__main__':
    main()