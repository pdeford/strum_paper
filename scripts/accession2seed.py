#!/usr/bin/env python

import hashlib
import struct
import sys

accession = sys.argv[1]
hashed = hashlib.sha1(accession).digest()
full_str = '{:0.8e}'.format(struct.unpack('d', hashed[:8])[0])
strip_str = full_str.replace('.','').replace('-','')

print strip_str[:8]