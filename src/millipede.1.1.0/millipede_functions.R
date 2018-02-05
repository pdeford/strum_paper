# MILLIPEDE functions
# Written by Kaixuan Luo. Duke University. 2012
# Citation: Luo, K. and Hartemink, A. J. Using DNase digestion data to accurately identify transcription factor binding sites. Pacific Symposium on Biocomputing, 2013.

# MILLIPEDE supervised version
# Parameters:
## data_train:        DNase I digestion data matrix to train the MILLIPEDE model for the transcription factor (training DNase data). The rows represent the PWM matched candidate binding sites, columns contain binned DNase cuts data. 
## data_test:         DNase I digestion data matrix for binding prediction using the trained model (testing DNase data). The rows represent the PWM matched candidate binding sites, columns contain binned DNase cuts data. 
## specificity_train: TF binding specificity information to train the MILLIPEDE model. e.g. training PWM scores.
# specificity_test:  TF binding specificity information for binding prediction using the trained model. e.g. testing PWM scores.
## label_train:       ChIP label to train the MILLIPEDE model. Positive sites have labels equal to one, negative sites have labels equal to zero.
# Output:
## millipede: a list of MILLIPEDE results containing
### prob: the binding probability predicted by MILLIPEDE model 
### model: the logistic regression model trained by MILLIPEDE

millipede = function(data_train, data_test, specificity_train, specificity_test, label_train){
  
  cat("Fit MILLIPEDE... \n")
  # Train MILLIPEDE model with training data and training labels
  data_train.df = data.frame(label = label_train, specificity = specificity_train, DNase = log2(data_train+0.1))
  model <- glm(formula = label ~ ., data = data_train.df, family = "binomial")
  # Predict the binding probability of the testing data using the model
  data_test.df = data.frame(specificity = specificity_test, DNase = log2(data_test+0.1))
  predictions = predict(model, data_test.df, type = "response", se = T)
  predictions.l = list(prob = predictions$fit, model = model)
  
  return(predictions.l)
  
}

# MILLIPEDE unsupervised version
# Parameters:
## data:        DNase I digestion data matrix for the transcription factor. The rows represent the PWM matched candidate binding sites, columns contain binned DNase cuts data. 
## specificity: TF binding specificity information. e.g. PWM score.
## weights_data:       provided weights for each bin of the data
## weights_specificity: provided weigths for TF binding specificity information. e.g. PWM score. Default weight for PWM score is 1.
## weights_intercept: provided intercept value. Default is 0, i.e. no intercept.
# Output:
## millipede_unsupervised: binding score (not probability) using the weights provided.

millipede_unsupervised = function(data, specificity, weights_data, weights_specificity = 1, weights_intercept = 0){
  
  cat("Fit MILLIPEDE unsupervised... \n")
  # Use the provided weights as the coefficents for covariates to get the prediction results
  data.matrix = as.matrix(cbind(1, specificity, DNase = log2(data+0.1)))  
  predict.logit = as.numeric(data.matrix %*% c(weights_intercept, weights_specificity, weights_data))
  predictions = 1 / (1 + exp(-predict.logit))
  
  return(predictions)
  
}

# MILLIPEDE partially unsupervised version
# Parameters:
## data:        DNase I digestion data matrix for the transcription factor. The rows represent the PWM matched candidate binding sites, columns contain binned DNase cuts data. 
## specificity: TF binding specificity information. e.g. PWM score.
## model:       MILLIPEDE model trained by this or other TFs
# Output:
## millipede_partial: binding probability predicted by MILLIPEDE model trained by this or other TFs.

millipede_partial = function(data, specificity, model){
  
  cat("Fit MILLIPEDE partially unsupervised... \n")
  # Use the reference model to predict the binding probability for the sites
  data.df = data.frame(specificity, DNase = log2(data+0.1))
  predictions = predict(model, data.df, type = "response", se = T)
  
  return(predictions$fit)
  
}

# Function for binning DNase cuts data with various binning methods
# Parameters:
## data:     DNase I digestion data matrix for the transcription factor. The rows represent the PWM matched candidate binding sites, columns contain 100bp flanking region on the left of motif, followed by motif region, and 100bp flanking on the right of motif. 
## bin_type: binning type for DNase cuts. Options for binning contain "M12", "M11", "M5", "M3", "M2", and "M1" as decribed in the paper. 
## FUN:      the way to represent the number of DNase cuts in each bin. Default option uses sum, other options such as mean or max could be specified.
# Output:
# bin_data: binned DNase cuts matrix

bin_data = function(data, bin_type, FUN = sum){
  
  motif_L = 101: round(ncol(data)/2)
  motif_R = (round(ncol(data)/2)+1) : (ncol(data)-100)
  
  flank_L1 = 81:100
  flank_L2 = flank_L1 - 20
  flank_L3 = flank_L2 - 20
  flank_L4 = flank_L3 - 20
  flank_L5 = flank_L4 - 20
  
  flank_R1 = (ncol(data)-99): (ncol(data)-80)
  flank_R2 = flank_R1 + 20
  flank_R3 = flank_R2 + 20
  flank_R4 = flank_R3 + 20
  flank_R5 = flank_R4 + 20
  
  M12 = data.frame(
    flank_L5 = apply(data[,flank_L5],1,FUN),
    flank_L4 = apply(data[,flank_L4],1,FUN),
    flank_L3 = apply(data[,flank_L3],1,FUN),
    flank_L2 = apply(data[,flank_L2],1,FUN),
    flank_L1 = apply(data[,flank_L1],1,FUN), 
    motif_L = apply(data[,motif_L],1,FUN),
    motif_R = apply(data[,motif_R],1,FUN),
    flank_R1 = apply(data[,flank_R1],1,FUN),
    flank_R2 = apply(data[,flank_R2],1,FUN),
    flank_R3 = apply(data[,flank_R3],1,FUN),
    flank_R4 = apply(data[,flank_R4],1,FUN),
    flank_R5 = apply(data[,flank_R5],1,FUN))
  
  left2_num = apply(M12[,c(1,2)],1,FUN)
  left1_num = apply(M12[,c(3,4,5)],1,FUN)
  motif_num = apply(M12[,c(6,7)],1,FUN)
  right1_num = apply(M12[,c(8,9,10)],1,FUN)
  right2_num = apply(M12[,c(11,12)],1,FUN)
  
  M11 = data.frame(M12[,c(1:5)], motif_num, M12[,c(8:12)])
  M5 = data.frame(left2_num, left1_num, motif_num, right1_num, right2_num)
  M3 = data.frame(left1_num, motif_num, right1_num)
  M2 = data.frame(left1_num, right1_num)
  M1 = left1_num + right1_num
  
  data_bin = switch(bin_type, 
                    M12 = M12, 
                    M11 = M11, 
                    M5 = M5,
                    M3 = M3,
                    M2 = M2,
                    M1 = M1)
  return(data_bin)
  
}
