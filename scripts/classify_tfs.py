#!/usr/bin/env python

import rdflib
from rdflib import Graph, URIRef, Literal


import regex as re

def remove_punctuation(text):
    return re.sub(ur"\p{P}+", "", text)
##################################################
# Create Mapping 
##################################################

g = Graph()
g.parse("data/tfclass.ttl", format="ttl")

def get_parents(s):
    parent_s = g.objects(s, subclass)
    for s2 in parent_s:
        if s2.toPython().startswith('http://sybig.de/tfclass'):
            label = g.label(s2)
            if label:
                label = label.toPython()
                break
    else:
        return []
    try:
        return [label] + get_parents(s2)
    except StopIteration:
        return [label]

genus = rdflib.term.Literal(u'Genus')
subclass = URIRef("http://www.w3.org/2000/01/rdf-schema#subClassOf")
superclass = URIRef("http://www.w3.org/2000/01/rdf-schema#superClassOf")
somevalues = URIRef("http://www.w3.org/2002/07/owl#someValuesFrom")
classification = {}
for s,p,o in g.triples( ( None, None, genus) ):
    tf = g.label(s).toPython().upper()
    # Check if it is found in humans
    for s2, p2, o2 in g.triples((s, None, None)):
        # if any(["Homo_sapiens" in z.toPython() for x,y,z in g.triples((o2, somevalues, None))]):
        classification[tf] = get_parents(s)[::-1]

##################################################
# Load Data
##################################################

f = open("output/correlations.txt")
tfs = []
for line in f:
    fields = line.split()
    tf = fields[0]
    tf = tf.split('.')[0]
    tf = tf.split('-')[-1]
    tfs.append(tf)

l = max([len(x) for x in tfs])

##################################################
# Extract Data
##################################################

classes = []
class_counts = {}

matches = []
loners = []

# Exact matches
for TF in tfs:
    tf = remove_punctuation(TF.upper())
    for key in classification.keys():
        if tf == remove_punctuation(key):
            matches.append((TF, key))
            break
    else:
        loners.append(TF)

# Partial matches, still high confidence
topop = []
for i, TF in enumerate(loners):
    tf = remove_punctuation(TF.upper())
    for key in classification.keys():
        if any([tf == bit for bit in remove_punctuation(key).split()]):
            matches.append((TF, key))
            topop.append(i)
            break
for i in topop[::-1]:
    loners.pop(i)

##################################################
# Store data for later use
##################################################

for tf, key in matches:
    nl = u"{}\t{}"
    pt2 = u"\t".join(classification[key][:4])
    print nl.format(tf, pt2).encode('utf-8')