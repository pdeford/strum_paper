#!/usr/bin/env python

import matplotlib as mpl
mpl.use("Agg")
import sys

try:
	import cPickle as pickle
except:
	import pickle


import learn_motifs
import logit

basename = sys.argv[1]
n_process = int(sys.argv[2])
random_seed = int(sys.argv[3])

models, sequences = learn_motifs.learn(basename, n_process, random_seed)
shuf_clf, flank_clf, shuf_results, flank_results = logit.score_all(basename, n_process, random_seed, models, sequences)

pickle.dump(list(models) + [shuf_clf, flank_clf,], open("output/{}.p".format(basename), 'wb'))
