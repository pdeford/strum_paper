#!/usr/bin/env bash

tf=$1
n=$(ls "data/${tf}."*".bed" |  cut -f2 -d '.' | sort | uniq | wc -l)

if [ $n == '1' ]; then
	# rm "data/${tf}."[^K]*
	echo "No other cell types for " $tf ". Exiting..."
	exit 5
fi

pattern="data/${tf}.K562.*.bed"
files=( $pattern )

# Record Accession numbers being used
echo "${files[0]}" >> output/dnase_accessions.txt
ls "data/${tf}."[^K]*".bed" >> output/dnase_accessions.txt

# Extract the peaks that are not in any of the other files
bedtools intersect -v -a "{files[0]}" -b "data/${tf}."[^K]*".bed" > "output/unique_${tf}.bed"

# Extract the peaks that are in every file
bedtools intersect -c -wa -a "{files[0]}" -b "data/${tf}."*".bed" | grep $n"$" > "output/ubiq_${tf}.bed"

# Store the peaks that are in neither all nor none of the other files
bedtools intersect -v -a "{files[0]}" -b "output/ubiq_${tf}.bed" "output/unique_${tf}.bed" > "output/med_${tf}.bed"

# Get some peaks that are not in K562
bedtools intersect -v -a $(ls "data/${tf}."[^K]*".bed" | head -n1) -b "{files[0]}" > "output/not_K562_${tf}.bed"