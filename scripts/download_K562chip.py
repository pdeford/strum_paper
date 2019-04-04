#!/usr/bin/env python

"""Programatically access and download all of the ChIP-seq data sets
from the human ENCODE project, subject to the following conditions:
- celltype : K562
- assay : ChIP-seq
- target : transcription factor
- assembly : hg19
- format : narrowpeak
- conservative IDR thresholded peaks

For each file associated with each transcription factor, rename to
comply with the convention: data/<factor>.K562.<accession>.bed
"""

import requests, json

def query(url):
	response = requests.get(url, headers=HEADERS)
	response_json_dict = response.json()
	return response_json_dict

def download_file(url, local_filename):
    r = requests.get(url, stream=True)
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024): 
            if chunk:
                f.write(chunk)
    return local_filename

celltype = "K562"

# Force return from the server in JSON format
HEADERS = {'accept': 'application/json'}


# This searches the ENCODE database
base = "https://www.encodeproject.org"
URL = base + "/search/?assay_title=ChIP-seq&target.investigated_as=transcription+factor&limit=all&assembly=hg19&biosample_ontology.term_name={}".format(celltype)

# Extract the JSON response as a python dict
response_json_dict = query(URL)

ids = []
targets = []
accessions = []
for thing in response_json_dict['@graph']:
	ids.append(thing['@id'])
	targets.append(thing['target']['label'])
	accessions.append(thing['accession'])

for i, target in enumerate(targets):
	response_json_dict = query(base + ids[i])
	for f in response_json_dict["files"]:
		if "assembly" in f:
			if "hg19" in f["assembly"] and "File" in f["@type"]:
				if "submitted_file_name" not in f: continue
			  	name = f["submitted_file_name"]
				if "narrowpeak.gz" in name.lower() and "conservative" in name.lower():
					href = f["href"]
					accession = href.split("/")[2]
					newname = ".".join([target, celltype, accession, "bed", "gz"])
					download_file(base + href, "data/" + newname)
					print "\t".join([target, celltype, accession])
					break
