#!/usr/bin/env python

import sys

try:
	import cPickle as pickle
except:
	import pickle

def main(tf):
	pwm, dwm, ml_strum, em_strum, (logit, scaler), (logit2, scaler2) = pickle.load(open("output/{}.p".format(tf), "rb"))
	coef = logit.coef_[0]

	print "\t".join([tf] + ["{:0.3f}".format(x) for x in coef])

if __name__ == '__main__':
	tf = sys.argv[1]
	main(tf)

