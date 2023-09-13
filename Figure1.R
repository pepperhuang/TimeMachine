library(limma)
library(glmnet)
source("TimeStampFns.R")

# This file contains expression data for all three datasets,
# time, subject, and condition metadata
load("data/TMexampleData.Rdata")

#=====================================================================
# convenience functions 
#=====================================================================

Rint <- function(...){
  Reduce(intersect,list(...))
}

asNum <- function(z){
  as.numeric(as.character(z))
}

nwhich <- function(x){
  if(is.null(names(x))){
    which(x)
  }else{
    names(which(x))
  }
}

asTime <- function(x){
  hrs <- (x%%24)%/%1; 
  min <- round(60*(x%%1)); 
  hrs[min==60] <- (hrs+1)%%24
  min[min==60] <- 0
  sprintf("%2i:%02i",hrs,min)
}

asTime0 <- function(x){
  hrs <- (x%%24)%/%1; 
  min <- round(60*(x%%1)); 
  hrs[min==60] <- (hrs+1)%%24
  min[min==60] <- 0
  sprintf("%02i:%02i",hrs,min)
}

#=====================================================================
# new plotting routines
#=====================================================================

timeErrPlot <- function(trueTimes,predTimes,...){
  errplot(trueTimes,predTimes,...)
  timeDiffs <- timeErr(trueTimes,predTimes)
  return(timeDiffs)
}

timeErrSignedPlot <- function(trueTimes,predTimes,...){
  errplot(trueTimes,predTimes,...)
  timeDiffs <- timeErrSigned(predTimes,trueTimes)
  return(timeDiffs)
}

tolplot <- function(trueHr,predHr,main="Time Since DLMO (24h)",add=FALSE,col=1,...){
  # plot % correct by tolerance, without lines 
  predTimes <- predHr%%24
  trueTimes <- trueHr%%24
  hrerr <- abs(predTimes-trueTimes)
  hrerr <- hrerr[!is.na(hrerr)]
  hrerr <- pmin(hrerr,24-hrerr)
  hrsoff <- seq(0,12,length=49)
  fracacc <- sapply(hrsoff,function(hrtol){
    100*sum(abs(hrerr)>hrtol)/length(hrerr)
  })
  if(!add){
    col=1
    plot(hrsoff,100-fracacc,xlim=c(0,12), type="n",
         main = main, xlab="",ylab="")
    mtext("correct to within (hrs)",side=1,line=2.2,cex=0.8)
    mtext(paste("% correct (N = ",length(hrerr),")",sep=""),side=2,line=2.4,cex=0.8)
    abline(a=0,b=100/12,col="grey")
    asTime <- function(x){
      hrs <- (x%%24)%/%1; 
      min <- round(60*(x%%1)); 
      if(min==60){hrs <- (hrs+1)%%24; min<-0}
      sprintf("%2i:%02i",hrs,min)
    }
  }
  lines(hrsoff, 100-fracacc, col=col,lwd=1.5,...)
  norm.fracacc <- (100-fracacc)/100
  norm.hrsoff <- hrsoff/12
  auc <- sum(norm.fracacc[-1]*diff(norm.hrsoff))
  return(list(auc=auc,mederr=median(abs(hrerr))) )
}

predplot <- function(trueHr,predHr,col=1,pch=1,main="Time of Day (24h)",...){
  # do both of the above -- set par(mfrow=c(2,1)) or (1,2) first!
  opar <- par(xpd=F,mar=c(4,4,3,1))	
  on.exit(par(opar))
  out <- timeErrPlot(trueHr,predHr,col,pch,main) 
  out <- tolplot(trueHr,predHr)
  invisible(out)
}

#=====================================================================
# For reproducibility
#=====================================================================
set.seed(194)
train.foldid <- sample(rep(seq(10),length=sum(all.meta$train)))

#=====================================================================
# Normalize Expression 
#=====================================================================
all.expr[is.na(all.expr)] <- 0
load("cycling_list.RData")
cycling.expr <- all.expr[rownames(all.expr) %in% cycling.list,]

## Compute ratio of genes or difference of log2 value
ratio.expr <- apply(cycling.expr,2,function(xx){
  tmp <- outer(xx,xx,'-')
  tmp[upper.tri(tmp)]
})

ratio.label <- apply(cycling.expr,2,function(xx){
  tmp <- outer(rownames(cycling.expr),rownames(cycling.expr),"paste")
  tmp[upper.tri(tmp)]
})

rownames(ratio.expr) <- ratio.label[,1]

## Z-Score
sample.mean <- apply(cycling.expr,2,mean)
sample.std <- apply(cycling.expr,2,sd)
predictor.expr.nor <- apply(cycling.expr, 1, function(x) {
  tmp <- (x-sample.mean)/sample.std
})
predictor.expr.nor <- t(predictor.expr.nor)
predictor.expr.final <- predictor.expr.nor

# select training samples 
CPtrain <- all.meta$train==1  & !is.na(all.meta$CPhrs)
all.foldid <- rep(0,length(all.meta$train))
all.foldid[all.meta$train==1] <- train.foldid
CPtrain.foldid <- all.foldid[CPtrain]
CPtrainRatioDat <- ratio.expr[,CPtrain]
CPtrainZScoreDat <- predictor.expr.final[,CPtrain]
CPtrainSubjs <- all.meta[CPtrain,"ID"]
CPtrainTimes <- all.meta[CPtrain,"CPhrs"]

# Train TimeMachine on Ratio
TM.CPhrs.Ratio <- trainTimeStamp(
  expr=CPtrainRatioDat, 
  subjIDs=CPtrainSubjs,
  times=CPtrainTimes,
  trainFrac=1, 
  recalib=FALSE,
  a = 0.2, s = 0.31,
  foldid=CPtrain.foldid, 
  plot=FALSE 
)
all.meta$TMpredRatio <- predTimeStamp(TM.CPhrs.Ratio,newx=ratio.expr, s = 0.31) 

# Train TimeMachine on ZScore
TM.CPhrs.ZScore <- trainTimeStamp(
  expr=CPtrainZScoreDat, 
  subjIDs=CPtrainSubjs,
  times=CPtrainTimes,
  trainFrac=1, 
  recalib=FALSE,
  a = 0.03, s = 0.3,
  foldid=CPtrain.foldid, 
  plot=FALSE 
)
all.meta$TMpredZScore <- predTimeStamp(TM.CPhrs.ZScore,newx=predictor.expr.final, s = 0.3) 

#=====================================================================
# PLOT Figure 
#=====================================================================
statLegend <- function(){
  legend("bottomright", bty="n", text.font=c(3,2,2,2), lwd=1.5, cex=1, lty=c(0,1,2,2), col=c("white","black","purple3","green4"), text.col=c("black","black","purple3","green4"),
         legend=c("         AUC,  med",sprintf(
           "%s: %.2f, %s",
           c(" rTM", " zTM", "   PL"),
           c(TMauc$auc, TMZauc$auc, PLSauc$auc),
           asTime(c(TMauc$mederr, TMZauc$mederr, PLSauc$mederr))
         )))
}


matPredDeg <- read.csv("matlab_PLSR_files/Prediction_PLSR_1sample_5_100_CPhrs.csv",header=F)
colnames(matPredDeg) <- c("CPdeg","matPLpred")
all.meta$PLmatpred <- matPredDeg$matPLpred*24/360

# pdf(file="Figure1.pdf",width=8,height=5.8)
par(mfcol=c(2,3))
# split up the results by study
all.meta.split <- split(all.meta, all.meta$study)

# Train/Test data (Moller) - TEST SAMPLES ONLY
TMauc <- with(all.meta.split$TrTe[all.meta.split$TrTe$train==0,],
              predplot(CPhrs,TMpredRatio,plot=T,main="Test Set (Moller et al.)",col=cond))
TMZauc <- with(all.meta.split$TrTe[all.meta.split$TrTe$train==0,],
               tolplot(CPhrs, TMpredZScore, add=T,col="purple3",lty=2))
PLSauc <- with(all.meta.split$TrTe[all.meta.split$TrTe$train==0,],
               tolplot(CPhrs, PLmatpred, add=T,col="green4",lty=2))
statLegend()

# V1 validation data (Archer)
TMauc <- with(all.meta.split$V1,
              predplot(CPhrs,TMpredRatio,plot=T, main="Validation V1 (Archer et al.)",col=cond))
TMZauc <- with(all.meta.split$V1,
               tolplot(CPhrs, TMpredZScore, add=T,col="purple3",lty=2))
PLSauc <- with(all.meta.split$V1,
               tolplot(CPhrs, PLmatpred, add=T,col="green4",lty=2))
statLegend()

# V2 validation data (Arnardottir &al)
# NO DLMO!

# V3 validation data (new RNA-seq)
TMauc <- with(all.meta.split$V3,
              predplot(CPhrs,TMpredRatio,plot=T,main="Validation V3 (RNA-seq)"))
TMZauc <- with(all.meta.split$V3,
               tolplot(CPhrs, TMpredZScore, add=T,col="purple3",lty=2))
PLSauc <- with(all.meta.split$V3,
               tolplot(CPhrs, PLmatpred, add=T,col="green4",lty=2))
statLegend()
# dev.off()
