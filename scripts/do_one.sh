#!/usr/bin/env bash

# Check command line arguments
if [[ $# -ne 2 ]]; then
    echo "ERROR: Illegal number of parameters:" $#
    echo "Usage:"
    echo "    ${0} n_process line"
    exit 5
fi

# Read command line arguments
n_process=$1
line=$2

# Extract TF name and accession
fname=${line##*/}
tf=${fname%%.*}
suffix=${fname##*.K562.}
accession=${suffix%.*}

basename=$tf'.'$accession
echo '---------------------------------------------------------------------'
echo -n "["$(date +"%F %T")"] "; echo $basename

# Set random seed for reproducibility
seed=$(python scripts/accession2seed.py $accession)
echo -n "["$(date +"%F %T")"] "; echo "Random seed:" $seed

# Sort on the signal value, from largest to smallest
echo -n "["$(date +"%F %T")"] "; echo "Sorting input"
sort -k7,7nr $line > $line'.sorted'

# Extract 100bp around the summit
echo -n "["$(date +"%F %T")"] "; echo "Find peak centers"
awk 'BEGIN{OFS="\t"}{a=($10+$2); print $1, a-50, a+50}' $line'.sorted' > 'data/'$basename'centered.bed'

# Prepare sequence data: Extract peaks and flanking sequences, run MEME
echo -n "["$(date +"%F %T")"] "; echo "Extract sequence for ChIP peaks"
bedtools getfasta -fi data/hg19.fa -bed 'data/'$basename'centered.bed' -fo 'data/'$basename'.fa'

echo -n "["$(date +"%F %T")"] "; echo "Get flanking  sequences as negative control"
bedtools flank -l 100 -r 100 -i $line'.sorted' -g data/hg19sizes.genome | head -n1000 > $line'.negative'
bedtools getfasta -fi data/hg19.fa -bed $line'.negative' -fo 'data/'$basename'.flank.fa'

echo -n "["$(date +"%F %T")"] "; echo "Run MEME"
meme-chip -oc 'output/'$basename'_meme' -dna -nmeme 500 -seed $seed -noecho \
  -norand -meme-nmotifs 1 -dreme-m 0 -spamo-skip \
  'data/'$basename'.fa'

# Learn each of the motifs: PWM, DWM, ML-StruM, EM-StruM
echo -n "["$(date +"%F %T")"] "; echo "Learn motifs and score seqs"
python scripts/do_peak_nonpeak.py $basename $n_process $seed >> output/chip_auc.txt

echo -n "["$(date +"%F %T")"] "; echo "Report coefficients of logit model"
./scripts/coeff.py $basename >> output/coefficents.txt

echo -n "["$(date +"%F %T")"] "; echo "Compare positioning of most significant matches"
./scripts/compare_positions.py $basename $n_process >> output/position_comp.txt

echo -n "["$(date +"%F %T")"] "; echo "Get correlation between motif scores"
./scripts/correlation.py $basename $seed >> output/correlations.txt