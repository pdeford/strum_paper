#!/usr/bin/env bash

# To track output and stderr:
##./do_all.sh 80 > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)

# Check command line arguments
if [[ $# -ne 1 ]]; then
    echo "ERROR: Illegal number of parameters:" $#
    echo "Usage:"
    echo "    ${0} n_process"
    exit 5
fi

n_process=$1

# Prepare the environment
source activate strum_paper
mkdir data figures output
touch output/correlations.txt

###########################
# Download necessary data #
###########################

## Download reference genome
#rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz ./data
#tar xvzf data/chromFa.tar.gz -C ./data
#cat data/chr*.fa > data/hg19.fa
#rm data/chr*.fa
#curl "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes" > data/hg19sizes.genome
#
## Download ChIP data for K562 cells
#python scripts/download_K562chip.py > data/accessions.txt
#
#gunzip data/*gz
#
## Download FOXA1 binding information from JASPAR
## Example for figures
#curl http://jaspar.genereg.net/download/sites/MA0148.1.sites > data/MA0148.1.sites

#############################################################
# Analyze each ChIP dataset for peak v non peak performance #
# and compare the performance of PWM to StruM               #
#############################################################

# Loop over each of the ChIP files
ls data/*K562*bed | while read line;
do
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

	bash scripts/do_one.sh $n_process $line
done

####################################
# Generate the figures from output #
####################################
python scripts/overview_fig.py
python scripts/generate_figures.py \
	-a output/chip_auc.txt \
	-r output/correlations.txt \
	-p output/position_comp.txt \
	-c output/coefficents.txt \
	> RESULTS.txt
