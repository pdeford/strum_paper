#!/usr/bin/env bash
:'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Execute this file with the first argument being !
! the number of cores available to be used.       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'

n_process=$1

# Prepare the environment
source activate strum_paper
mkdir data figures output
touch output/correlations.txt

###########################
# Download necessary data #
###########################

# Download reference genome
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz ./data
tar xvzf data/chromFa.tar.gz -C ./data
cat data/chr*.fa > data/hg19.fa
rm data/chr*.fa
curl "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes" > data/hg19sizes.genome

# Download ChIP data for K562 cells
python scripts/download_K562chip.py > data/accessions.txt

gunzip data/*gz

# Download FOXA1 binding information from JASPAR
# Example for figures
curl http://jaspar.genereg.net/download/sites/MA0148.1.sites > data/MA0148.1.sites

#############################################################
# Analyze each ChIP dataset for peak v non peak performance #
# and compare the performance of PWM to StruM               #
#############################################################

# Set a random seed for reproducibility
declare -i seed=602

# Loop over each of the ChIP files
ls data/*K562*bed | while read line;
do
	# Iterate seed
	seed+=1

	# Extract TF name and accession
	fname=${line##*/}
	tf=${fname%%.*}
	suffix=${fname##*.K562.}
	accession=${suffix%.*}

	basename=$tf'.'$accession

	# Check if this iteration has already been done. If so,
	##continue to the next one. (Allows several versions to
	##run simultaneously.)
	if [ -e 'data/'$basename'.fa' ]; then
		continue
	fi

	echo -n "["$(date +"%F %T")"] "; echo $basename

	# Prepare sequence data: Extract peaks and flanking sequences, run MEME
	echo -n "["$(date +"%F %T")"] "; echo "Extract sequence for ChIP peaks"
	bedtools getfasta -fi data/hg19.fa -bed $line -fo 'data/'$basename'.fa'

	echo "Get flanking  sequences as negative control"
	bedtools flank -l 0.0 -r 1.0 -pct -i $line -g data/hg19sizes.genome > $line'.negative'
	bedtools getfasta -fi data/hg19.fa -bed $line'.negative' -fo 'data/'$basename'.flank.fa'

	echo -n "["$(date +"%F %T")"] "; echo "Run MEME"
	meme-chip -oc 'output/'$basename'_meme' -dna -nmeme 500 -seed $seed -noecho \
	  -norand -meme-nmotifs 1 -dreme-m 1 -spamo-skip 'data/'$basename'.fa'

	# Learn each of the motifs: PWM, DWM, ML-StruM, EM-StruM
	echo -n "["$(date +"%F %T")"] "; echo "Learn motifs and score seqs"
	python scripts/do_peak_nonpeak.py $basename $n_process $seed >> output/chip_auc.txt

	echo -n "["$(date +"%F %T")"] "; echo "Report coefficients of logit model"
	./scripts/coeff.py $basename >> output/coefficents.txt

	echo -n "["$(date +"%F %T")"] "; echo "Compare positioning of most significant matches"
	./scripts/compare_positions.py $basename $n_process >> output/position_comp.txt

done

##############################################################
# Find the correlation between scores assigned by each model #
##############################################################

# Loop over each of the ChIP files
ls data/*K562*bed | while read line;
do
	# Extract TF name and accession
	fname=${line##*/}
	tf=${fname%%.*}
	suffix=${fname##*.K562.}
	accession=${suffix%.*}

	basename=$tf'.'$accession

	echo "Get correlation between motif scores"
	./scripts/correlation.py $basename >> output/correlations.txt
done;

####################################
# Generate the figures from output #
####################################
python overview_fig.py
python scripts/generate_figures.py \
	-a output/chip_auc.txt \
	-r output/correlations.txt \
	-p output/position_comp.txt \
	-c output/coefficents.txt
