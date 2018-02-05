#!/usr/bin/env bash

tf=$1
n=$(ls "data/${tf}."*".bed" | wc -l)

if ls "data/${tf}.K562."*".bed" 1> /dev/null 2>&1; then
	rm "data/${tf}"*
	exit 5
else
	if [ n == '1' ]; then
		rm "data/${tf}."[^K]*
		exit 5
	fi
fi

# Merge all K562 files for TF
cat "data/${tf}.K562."*".bed" | sortBed -i stdin | mergeBed -i stdin > "data/${tf}.K562.merged.bed"

# Extract the peaks that are not in any of the other files
bedtools intersect -v -a "data/${tf}.K562.merged.bed" -b "data/${tf}."[^K]*".bed" > "output/unique_${tf}.bed"

# Extract the peaks that are in every file
bedtools intersect -c -wa -a "data/${tf}.K562.merged.bed" -b "data/${tf}."*".bed" | grep $n"$" > "output/ubiq_${tf}.bed"

# Store the peaks that are in neither all nor none of the other files
bedtools intersect -v -a "data/${tf}.K562.merged.bed" -b "output/ubiq_${tf}.bed" "output/unique_${tf}.bed" > "output/med_${tf}.bed"

# Get some peaks that are not in K562
bedtools intersect -v -a $(ls "data/${tf}."[^K]*".bed" | head -n1) -b "data/${tf}.K562.merged.bed" > "output/not_K562_${tf}.bed"

