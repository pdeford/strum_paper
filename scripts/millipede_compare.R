args <- commandArgs(TRUE)
factor <- args[1]

source("millipede.1.1.0/millipede_functions.R")
library(ROCR)
library(pracma)

# LOAD PERFORMANCE OF STRUM MODEL FROM `strum_dnase.py`
auc1 <- read.table(paste("output/auc1_", factor, ".txt", sep=""), header=F, strip.white=T)
auc1 <- data.matrix(auc1)

auc2 <- read.table(paste("output/auc2_", factor, ".txt", sep=""), header=F, strip.white=T)
auc2 <- data.matrix(auc2)

auc3 <- read.table(paste("output/auc3_", factor, ".txt", sep=""), header=F, strip.white=T)
auc3 <- data.matrix(auc3)

# LOAD DATA OUTPUT FROM `strum_dnase.py`
train_labels <- read.table(paste("output/train_label_", factor, ".txt", sep=""), header=F, strip.white=T)
train_labels <- c(data.matrix(train_labels))

test_labels <- read.table(paste("output/test_label_", factor, ".txt", sep=""), header=F, strip.white=T)
test_labels <- c(data.matrix(test_labels))

train_motifs <- read.table(paste("output/train_motif_", factor, ".txt", sep=""), header=F, strip.white=T)
train_motifs <- c(data.matrix(train_motifs))

test_motifs <- read.table(paste("output/test_motif_", factor, ".txt", sep=""), header=F, strip.white=T)
test_motifs <- c(data.matrix(test_motifs))

train_dnase <- read.table(paste("output/train_dnase_", factor, ".txt", sep=""), header=F, strip.white=T)
train_dnase <- data.matrix(train_dnase)

test_dnase <- read.table(paste("output/test_dnase_", factor, ".txt", sep=""), header=F, strip.white=T)
test_dnase <- data.matrix(test_dnase)

############################################################

train_pwm <- read.table(paste("output/train_pwm_stuff_", factor, ".txt", sep=""), header=F, strip.white=T)
train_pwm <- c(data.matrix(train_pwm))

test_pwm <- read.table(paste("output/test_pwm_stuff_", factor, ".txt", sep=""), header=F, strip.white=T)
test_pwm <- c(data.matrix(test_pwm))

train_pwm_dnase <- read.table(paste("output/train_pwm_dnase_", factor, ".txt", sep=""), header=F, strip.white=T)
train_pwm_dnase <- data.matrix(train_pwm_dnase)

test_pwm_dnase <- read.table(paste("output/test_pwm_dnase_", factor, ".txt", sep=""), header=F, strip.white=T)
test_pwm_dnase <- data.matrix(test_pwm_dnase)

############################################################

# Bin DNase data
train_dnase <- bin_data(train_dnase, "M5") 
test_dnase <- bin_data(test_dnase, "M5")
train_pwm_dnase <- bin_data(train_pwm_dnase, "M5") 
test_pwm_dnase <- bin_data(test_pwm_dnase, "M5")

# Pull some negative examples to use for training the model
n1 <- length(train_labels)
n2 <- length(test_labels)
train_labels2 <- c(train_labels[1:n1], test_labels[(n2-n1+1):n2])
train_motifs2 <- c(train_motifs[1:n1], test_motifs[(n2-n1+1):n2])
train_pwm2 <- c(train_pwm[1:n1], test_pwm[(n2-n1+1):n2])
train_dnase2 <- rbind(train_dnase[1:n1,], test_dnase[(n2-n1+1):n2,])
train_pwm_dnase2 <- rbind(train_pwm_dnase[1:n1,], test_pwm_dnase[(n2-n1+1):n2,])

# Train model with MILLIPEDE
data_train <- train_dnase2
data_test <- test_dnase
specificity_train <- train_motifs2
specificity_test <- test_motifs
label_train <- train_labels2

millipede_result <- millipede(data_train, data_test, specificity_train, specificity_test, label_train) # supervised MILLIPEDE learning
pred_milli <- millipede_result$prob # predicted probability of binding from MILLIPEDE
model_milli <- millipede_result$model # model trained by MILLIPEDE

millipede_result2 <- millipede(train_pwm_dnase2, test_pwm_dnase, train_pwm2, test_pwm, label_train) # supervised MILLIPEDE learning
pred_milli2 <- millipede_result2$prob # predicted probability of binding from MILLIPEDE
model_milli2 <- millipede_result2$model # model trained by MILLIPEDE

# Check performance
pred <- prediction(pred_milli, test_labels)
perf <- performance(pred, "tpr", "fpr")
#plot(perf)

perf.auc <- performance(pred, 'auc')
auc <- perf.auc@y.values
print(auc)

pred2 <- prediction(pred_milli2, test_labels)
perf2 <- performance(pred2, "tpr", "fpr")

perf2.auc <- performance(pred2, 'auc')
auc2b <- perf2.auc@y.values
print(auc2b)

# Plot comparison to StruM performance

AUCa <- trapz(auc3[1,], auc3[2,])
AUCb <- trapz(auc1[1,], auc1[2,])
AUCe <- trapz(auc2[1,], auc2[2,])


legend.text <- c(
  paste("PWM", " (", formatC(AUCa[[1]], digits=2, format="f"), ")", sep=""),
  paste("PWM + Millipede", " (", formatC(auc2b[[1]], digits=2, format="f"), ")", sep=""),
  paste("StruM", " (", formatC(AUCb[[1]], digits=2, format="f"), ")", sep=""),
  paste("StruM + Millipede", " (", formatC(auc[[1]], digits=2, format="f"), ")", sep=""),
  paste("DNase-StruM", " (", formatC(AUCe[[1]], digits=2, format="f"), ")", sep="")
)

pdf(paste("output/mill_comp_", factor, ".pdf", sep=""))
plot(c(0,1), c(0,1), type='l', xlab=perf@x.name, ylab=perf@y.name, 
     ylim=c(0,1), xlim=c(0,1), col='grey', lty=2)
lines(auc3[1,], auc3[2,], ylim=c(0,1), xlim=c(0,1), lty='dotted', col='orange')
lines(perf2@x.values[[1]], perf2@y.values[[1]], ylim=c(0,1), xlim=c(0,1), lty='longdash', col='orange')
lines(auc1[1,], auc1[2,], ylim=c(0,1), xlim=c(0,1), lty='dotted', col='blue')
lines(perf@x.values[[1]], perf@y.values[[1]], ylim=c(0,1), xlim=c(0,1), lty='longdash', col='blue')
lines(auc2[1,], auc2[2,], ylim=c(0,1), xlim=c(0,1), lty="solid", col='blue')
legend("bottomright", legend=legend.text,
       col=c('orange', 'orange', 'blue', 'blue', 'blue' ), lty=c(3, 5, 3, 5, 1))
dev.off()


fileConn <- file(paste("output/", factor, "_AUCs.txt", sep=""))
writeLines(legend.text, fileConn)
close(fileConn)
