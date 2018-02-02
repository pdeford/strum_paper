#!/usr/bin/env python

import sys
for line in open(sys.argv[1]):
    fields = line.split()
    end = int(fields[2])
    l = end-int(fields[1])
    print "{}\t{}\t{}".format(fields[0], end, end+l)