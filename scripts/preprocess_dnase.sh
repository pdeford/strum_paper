#!/usr/bin/env bash

# Usage:
# ls data/*.bed | cut -d "." -f1 | cut -d "/" -f2 | sort | uniq | while read tf; do bash scripts/preprocess.sh "$tf"; done

tf=$1
n=$(ls "data/${tf}."*".bed" | wc -l)

if [ ! -e "data/${tf}.K562."* ]; then
	rm "data/${tf}"*
	exit 5
else
	if [ n == '1' ]; then
		rm "data/${tf}."[^K]*
		exit 5
	fi
fi


# Extract the peaks that are not in any of the other files
bedtools intersect -v -a "data/${tf}.K562."* -b "data/${tf}."[^K]* > "output/unique_${tf}.bed"

# Extract the peaks that are in every file
bedtools intersect -c -wa -a "data/${tf}.K562."* -b "data/${tf}."*".bed" | grep $n"$" > "output/ubiq_${tf}.bed"

# Store the peaks that are in neither all nor none of the other files
bedtools intersect -v -a "data/${tf}.K562."* -b "output/ubiq_${tf}.bed" "output/unique_${tf}.bed" > "output/med_${tf}.bed"

# Get some peaks that are not in K562
bedtools intersect -v -a $(ls "data/${tf}."[^K]* | head -n1) -b "data/${tf}.K562."* > "output/not_K562_${tf}.bed"