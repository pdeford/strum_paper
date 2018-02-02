#!/usr/bin/env bash
:'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Execute this file with the first argument being !
! the number of cores available to be used.       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'

n_process=$1

# Prepare the environment
source activate strum
mkdir data figures output

###########################
# Download necessary data #
###########################

# Download reference genome
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz ./data
tar xvzf data/chromFa.tar.gz -C ./data
cat data/chr*.fa > data/hg19.fa
#rm data/chr*.fa

# Download ChIP data for K562 cells
python scripts/download_K562chip.py > data/accessions.txt
# Download additional ChIP data for cell type specific predictions
python scripts/download_OTHERchip.py >> data/accessions.txt

gunzip data/*gz

# Download DNase data for K562 cells
dnase_accession="ENCFF111KJD"
curl -o "data/DNase.K562."$dnase_accession".bw" \
	"https://www.encodeproject.org/files/"$dnase_accession"/@@download/"$dnase_accession".bigWig"
echo $'DNASE-seq\tK562\t'$dnase_accession >> data/accessions.txt

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

	echo $basename

	# Prepare sequence data: Extract peaks and flanking sequences, run MEME
	echo "Extract sequence for ChIP peaks"
	bedtools getfasta -fi data/hg19.fa -bed $line -fo 'data/'$basename'.fa'

	echo "Get flanking  sequences as negative control"
	python -c get_flanks.py $line > $line'.negative'
	bedtools getfasta -fi data/hg19.fa -bed $line'.negative' -fo 'data/'$basename'.flank.fa'

	echo "Run MEME"
	meme-chip -oc 'output/'$basename'_meme' -dna -nmeme 500 -seed $seed -noecho \
	  -norand -meme-nmotifs 1 -dreme-m 1 -spamo-skip 'data/'$basename'.fa'

	# Learn each of the motifs: PWM, DWM, ML-StruM, EM-StruM
	python scripts/do_peak_nonpeak.py $basename $n_process $seed >> output/chip_auc.txt

	echo "Report coefficients of logit model"
	./scripts/coeff.py $tf >> output/coefficents.txt

	echo "Compare positioning of most significant matches"
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

####################################################################
# Analyze performance on cell type specific predicitons with DNase #
####################################################################

seedstart=711

function inside_loop {
	tf=$1
	echo $tf
	bash scripts/preprocess_dnase.sh $tf

	if [ $? != 5 ]; then
		# # EXTRACT SEQUENCES FROM BED TO FA
		bedtools getfasta -fo "output/${tf}.fa" -fi data/hg19.fa -bed "output/ubiq_${tf}.bed"
		# # RUN MEME ON EXTRACTED SEQUENCES
		meme-chip -oc "output/${tf}_meme" -dna -order 2 -nmeme 1000 -seed $seedstart -norand -noecho -meme-nmotifs 1 -dreme-m 1 "output/${tf}.fa"
		# RUN STRUM_DNASE + MILLIPEDE
		python scripts/strum_dnase.py "data/DNase.K562.ENCFF111KJD.bw" data/ $tf
		Rscript scripts/millipede_compare.R $tf
	fi
}

# Extract sequence, run MEME, do cell type specific analysis
ls data/*.bed | cut -d "." -f1 | cut -d "/" -f2 | sort | uniq |\
	xargs -n 1 -P $n_process inside_loop

# Consolidate the output
ls output/*_AUCs.txt | while read line;
do
	fname=${line##*/}
	name=${fname%%.*}
	tf=${name%%_*}
	echo -n $tf$'\t'
	awk '{printf "%s\t",substr($NF,2,4)}' 'output/'$tf'_AUCs.txt'
	echo
done > output/dnase_consolidated.txt

####################################
# Generate the figures from output #
####################################
python overview_fig.py
python scripts/generate_figures.py \
	output/chip_auc.txt \
	output/correlations.txt \
	output/dnase_consolidated.txt \
	output/position_comp.txt
