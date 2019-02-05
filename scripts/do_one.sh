
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
echo -n "["$(date +"%F %T")"] "; echo $basename

# Set random seed for reproducibility
seed=$(python scripts/accession2seed.py $accession)

# Extract 100bp around the summit
sort -k7,7nr $line | awk 'BEGIN{OFS="\t"}{a=($10+$2); print $1, a-50, a+50}' > 'data/'$basename'centered.bed'

# Prepare sequence data: Extract peaks and flanking sequences, run MEME
echo -n "["$(date +"%F %T")"] "; echo "Extract sequence for ChIP peaks"
bedtools getfasta -fi data/hg19.fa -bed 'data/'$basename'centered.bed' -fo 'data/'$basename'.fa'

echo "Get flanking  sequences as negative control"
bedtools flank -l 0.0 -r 1.0 -pct -i $line -g data/hg19sizes.genome > $line'.negative'
bedtools getfasta -fi data/hg19.fa -bed $line'.negative' -fo 'data/'$basename'.flank.fa'

echo -n "["$(date +"%F %T")"] "; echo "Run MEME"
meme-chip -oc 'output/'$basename'_meme' -dna -nmeme 500 -seed $seed -noecho \
  -norand -meme-nmotifs 1 -dreme-m 0 -spamo-skip 'data/'$basename'.fa'

# Learn each of the motifs: PWM, DWM, ML-StruM, EM-StruM
echo -n "["$(date +"%F %T")"] "; echo "Learn motifs and score seqs"
python scripts/do_peak_nonpeak.py $basename $n_process $seed >> output/chip_auc.txt

echo -n "["$(date +"%F %T")"] "; echo "Report coefficients of logit model"
./scripts/coeff.py $basename >> output/coefficents.txt

echo -n "["$(date +"%F %T")"] "; echo "Compare positioning of most significant matches"
./scripts/compare_positions.py $basename $n_process >> output/position_comp.txt

echo -n "["$(date +"%F %T")"] "; echo "Get correlation between motif scores"
./scripts/correlation.py $basename $seed >> output/correlations.txt