# MILLIPEDE example
# Written by Kaixuan Luo. Duke University. 2012
# Citation: Luo, K. and Hartemink, A. J. Using DNase digestion data to accurately identify transcription factor binding sites. Pacific Symposium on Biocomputing, 2013.

# Data preparation:
# DNaseData contains the DNase I digestion data for each transcription factor can be preapred as a matrix. The rows represent the PWM matched sites, columns contain 100bp flanking region on the left of motif, followed by motif region, and 100bp flanking on the right of motif. 
# pwm contains the PWM scores matching to the motif.
# label contains the ChIP labels from ChIP experiments. Positive sites have labels equal to one, negative sites have labels equal to zero.

# Source MILLIPEDE functions
source("millipede_functions.R")

# Load DNase digestion data, TF binding specificity and ChIP labels
load("Rap1.Rdata")

# Select training sites and test sites. 
# E.g. choose sites of 1 to 1000 as test sites, and the rest sites are training sites. 
testSite = 1:1000
trainSite = -testSite

data_train <- DNaseData[trainSite,]
data_test <- DNaseData[testSite,]
specificity_train <- pwm[trainSite]
specificity_test <- pwm[testSite]
label_train <- label[trainSite]
label_test <- label[testSite]

# Bin the DNase cuts data
# The bin_data function count the DNase cuts number within each bin. Options for binning contain "M12", "M11", "M5", "M3", "M2", and "M1" as decribed in the paper. 
data_train <- bin_data(data_train, "M5") 
data_test <- bin_data(data_test, "M5")

# MILLIPEDE supervised version
millipede_result <- millipede(data_train, data_test, specificity_train, specificity_test, label_train) # supervised MILLIPEDE learning
pred_milli <- millipede_result$prob # predicted probability of binding from MILLIPEDE
model_milli <- millipede_result$model # model trained by MILLIPEDE
# save(model_milli, file = "Rap1.M5.model")

# MILLIPEDE unsupervised version
weights_data = c(1, 1, -1, 1, 1) # set the weights for DNase bins in M5 model, other choice of weights could be chosen according to users' preference. 
weights_specificity = 1 # The default sets the weight for PWM score to 1, other choice of weight could be chosen according to users' preference.
pred_milli.u <- millipede_unsupervised(data_test, specificity_test, weights_data, weights_specificity) # unsupervised MILLIPEDE predictions

# MILLIPEDE partially unsupervised version
Reb1.M5.model <- get(load("Reb1.M5.model")) # Load the MILLIPEDE model trained by this or other TFs
pred_milli.p <- millipede_partial(data_test, specificity_test, Reb1.M5.model) # partially unsupervised MILLIPEDE predictions

cat("MILLIPEDE example finished.\n")