RNGversion("3.5.1")
library(limma)
library(glmnet)
source("TimeStampFns.R")

# This file contains expression data for all three datasets,
# time, subject, and condition metadata, as well as previously
# obtained ZeitZeiger and PLSR predictions for comparison
load("data/TMexampleData.Rdata")

#=====================================================================
# PLSR
#=====================================================================
# Normalized the data
all.expr[is.na(all.expr)] <- 0

library(matrixStats)
train.expr <- all.expr[,all.meta$train == 1]
train.meta <- all.meta[all.meta$train == 1,]

# Qunatile Normalized Based on Training Data
df_sorted <- data.frame(apply(train.expr,2,sort))
df_mean <- apply(df_sorted,1,mean)

index_to_mean <- function(my_index, my_mean){
  return(my_mean[my_index])
}

df_rank <- apply(all.expr,2,rank,ties.method = "min")
all.expr.final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
rownames(all.expr.final) <- rownames(all.expr)


library(R.matlab)
CPtrain <- all.meta$train==1 & !is.na(all.meta$CPhrs)
all.meta$LocalDeg <- 360*(all.meta$LocalTime/24)
# write out training data
time <- all.meta[CPtrain,"CPdeg"]
num <- t(all.expr.final[,CPtrain])
writeMat("matlab_PLSR_files/CPtrain_CPhrs.mat",dat=list(num=num,time=time))
# write out validation data
time <- all.meta[,"CPdeg"]
num <- t(all.expr.final) # not standardized
writeMat("matlab_PLSR_files/CPall_CPhrs.mat",dat=list(num=num,time=time))


# write out training data
CPtrain <- all.meta$train==1 & !is.na(all.meta$LocalDeg)
time <- all.meta[CPtrain,"LocalDeg"]
num <- t(all.expr.final[,CPtrain])
writeMat("matlab_PLSR_files/CPtrain_LocalTime.mat",dat=list(num=num,time=time))
# write out validation data
time <- all.meta[,"LocalDeg"]
num <- t(all.expr.final) # not standardized
writeMat("matlab_PLSR_files/CPall_LocalTime.mat",dat=list(num=num,time=time))


