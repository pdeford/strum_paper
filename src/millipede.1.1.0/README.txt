MILLIPEDE is a model for identifying transcription factor binding sites. MILLIPEDE is described in the paper: Luo, K. and Hartemink, A. J. Using DNase digestion data to accurately identify transcription factor binding sites. Pacific Symposium on Biocomputing, 2013.

MILLIPEDE was developed by Kaixuan Luo under the direction of Alexander J. Hartemink in the Program in Computational Biology and Bioinformatics and Department of Computer Science at Duke University.

MILLIPEDE model contains supervised, unsupervised and partially unsupervised versions. The logistic regression framework permits flexible ways to include information from different experimental resources.

Supervised MILLIPEDE builds logistic regression models using binned DNase data and TF binding specificity information as covariates trained by ChIP labels. When ChIP labels for the particular TF in the particular condition are not available, unsupervised and partially unsupervised versions can be used as described in the paper. 
 
Data preparation:
Supervised MILLIPEDE requires binned DNase data, TF binding specificity information (mainly using PWM score), ChIP labels (0 or 1) as in the example data.

DNase data matrix could be binned using the bin_data function in millipede_functions.R, or in users' preferred ways. After binning, the DNase data should be repented as a matrix, in which rows represents the candidate binding sites, and columns represents the number of DNase cuts within each bin. 

MILLIPEDE was originally written in R. The current version of MILLIPEDE is 1.1.0.
millipede_functions.R contains the main functions for MILLIPEDE models and the bin_data function to bin DNase cuts data.

millipede_example.R shows an example to bin the DNase data, and run MILLIPEDE models in supervised, unsupervised and partially unsupervised versions.

For questions or more information, please contact kaixuan.luo@duke.edu

