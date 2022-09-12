
getInitParameters <- function(segmVect,meanVect,sdVect, truncBoundRight,truncBoundLeft){
  
  truncBoundLeftVect<-c(meanVect[2]*3,meanVect[2]-(3*sdVect[2]),-truncBoundLeft,truncBoundRight,meanVect[4]+(3*sdVect[4]))
  truncBoundRightVect<-c(meanVect[2]-(3*sdVect[2]),-truncBoundLeft,truncBoundRight,meanVect[4]+(3*sdVect[4]),meanVect[4]*3)
  
  for (i in 1:length(meanVect)){ 
    ind <- which(segmVect<=truncBoundRightVect[i] & segmVect>=truncBoundLeftVect[i])
    
    if (length(ind)==1){
      meanVect[i]<-segmVect[ind]
    }else if (length(ind)>1){
      meanVect[i]<-mean(segmVect[ind])
      sdVect[i]<-sd(segmVect[ind])
    }
  }
  
  truncBoundLeftVect<-c(meanVect[2]*3,meanVect[2]-(3*sdVect[2]),-truncBoundLeft,truncBoundRight,meanVect[4]+(3*sdVect[4]))
  truncBoundRightVect<-c(meanVect[2]-(3*sdVect[2]),-truncBoundLeft,truncBoundRight,meanVect[4]+(3*sdVect[4]),meanVect[4]*3)
  
  InitParameters<-list()
  InitParameters$meanVect<-meanVect
  InitParameters$sdVect<-sdVect
  InitParameters$truncBoundLeftVect<-truncBoundLeftVect
  InitParameters$truncBoundRightVect<-truncBoundRightVect
  
  InitParameters
}

####### Expectation Maximization (EM) algorithm #############
EM <- function(segmVect,meanVect,sdVect,truncBoundRight,truncBoundLeft, prior){
  
  InitParameters<-getInitParameters(segmVect,meanVect,sdVect,truncBoundRight,truncBoundLeft)
  meanVect<-InitParameters$meanVect
  sdVect<-InitParameters$sdVect
  truncBoundLeftVect<-InitParameters$truncBoundLeftVect
  truncBoundRightVect<-InitParameters$truncBoundRightVect
  
  LikeliNew<-sum(getPosteriorProbability(segmVect,meanVect,sdVect,prior)*prior)
  
  threshold<-1e-03
  for (i in 1:100){
    meanVectold<-meanVect
    LikeliOld<-LikeliNew
    taux<-EStep(segmVect,meanVect,sdVect,prior,truncBoundLeftVect,truncBoundRightVect)
    MResult<-MStep(segmVect,taux,meanVect,sdVect)
    meanVect<-MResult$mu
    sdVect<-MResult$sdev
    prior<-MResult$prior
    
    LikeliNew<-sum(getPosteriorProbability(segmVect,meanVect,sdVect,prior)*prior)
    if (abs(LikeliNew-LikeliOld)<threshold){	
      break
    }
    
  }
  ResultsEM<-list()
  ResultsEM$meanVect<-meanVect
  ResultsEM$sdVect<-sdVect
  ResultsEM$prior<-prior
  ResultsEM$iter<-i
  ResultsEM$bound<-cbind(truncBoundLeftVect,truncBoundRightVect)
  ResultsEM
}

######## Truncated Gaussian Density  Function ########## EQ 16

getTruncatedGaussianDensity <- function(x,mean,sd,truncBoundLeft,truncBoundRight){
  
  c((dnorm(x,mean,sd)*(x<=truncBoundRight)*(x>=truncBoundLeft))/(pnorm(truncBoundRight, mean, sd)-pnorm(truncBoundLeft,mean,sd)))
}

######## Conditional Probabilities Function ########## EQ (18)
getConditionalProbabilities <-function(meanVect,prior,normaldataMat){
  
  do.call(cbind,lapply(1:length(meanVect), function(x) c(prior[x]*(t(normaldataMat[,x])))))
}

######## Posterior Probability Function ##########
getPosteriorProbability<-function(segmVect,meanVect,sdVect,prior){
  
  normaldataMat <- mapply(function(x,y) dnorm(segmVect,mean=x,sd=y), meanVect, sdVect)
  
  tauMat <- getConditionalProbabilities(meanVect,prior,normaldataMat)
  
  Posterior <- tauMat/rowSums(tauMat)
  Posterior
}

####### Expectation Step #############
EStep<-function(segmVect,meanVect,sdVect,prior,truncBoundLeftVect,truncBoundRightVect){
  
  normaldataMat <- do.call(cbind,lapply(1:length(meanVect), function(x) {
    mu<-meanVect[x]
    sdev<-sdVect[x]
    l<-truncBoundLeftVect[x]
    u<-truncBoundRightVect[x]
    if (sum((segmVect<=u)*(segmVect>=l))!=0){
      
      normaldataVec<- getTruncatedGaussianDensity(t(segmVect),mu,sdev,l,u) #EQ 16
      
      if(any(is.infinite(normaldataVec))){
        normaldataVec[which(is.infinite(normaldataVec))]<-100
      }
      
    }else{
      normaldataVec <- rep(0,length(segmVect))
    }
    normaldataVec
  }))
  
  
  tauMat <- getConditionalProbabilities(meanVect,prior,normaldataMat) #EQ (18)
  
  deno <- rowSums(tauMat)
  ind0 <- which(deno==0)
  
  if (length(ind0)!=0){
    for (ind in ind0){
      indmin <- which.min(abs(meanVect-segmVect[ind]))
      tauMat[ind,indmin] <- 1
      deno[ind] <- 1
    }
  }
  
  taux<- tauMat/deno
  taux
}

########### Maximization Step ###################
MStep<-function(segmVect,taux,meanVect,sdVect){
  
  mean <- meanVect
  sdev <- sdVect
  priorTemp <- colSums(taux)
  prior <- rep(1e-06, length(meanVect))
  for (j in 1:ncol(taux)){
    if (priorTemp[j]!=0){
      mean[j] <- (taux[,j]%*%segmVect)/priorTemp[j] #EQ 19
      sdevNew <- sqrt((taux[,j]%*%((segmVect-mean[j])^2))/priorTemp[j]) #EQ 20
      if (sdevNew > 0){ 
        sdev[j]<- sdevNew
      }
      prior[j] <- priorTemp[j]/length(segmVect) #EQ 21
    } 
  }
  
  prior <- prior/sum(prior)
  ParamResult <- list()
  ParamResult$mu <- mean
  ParamResult$sdev <- sdev
  ParamResult$prior <- prior
  ParamResult
}


getExtractSegm <- function(count_mtx, breaks, par_cores = 20){
  
  n <- nrow(count_mtx)
  
  startBreaks <- breaks[-length(breaks)]
  endBreaks <- c(startBreaks[-1])
  endBreaks <- c(endBreaks-1,n)
  
  funcCNA <- function(cell){
    
    x<-c()
    for (i in 1:length(startBreaks)){
      
      meanValue <- mean(count_mtx[startBreaks[i]:endBreaks[i],cell])
      
      x <- rbind(x,c(cell,startBreaks[i],endBreaks[i],meanValue,0))
      
    }
    return(x)
  }
  
  if(Sys.info()["sysname"]=="Windows"){
    
    cl <- parallel::makeCluster(getOption("cl.cores", par_cores))
    
    seg.test <- parallel::parLapply(cl, 1:ncol(count_mtx), funcCNA)
    
    parallel::stopCluster(cl)
  }else{
    seg.test <- parallel::mclapply(1:ncol(count_mtx), funcCNA, mc.cores = par_cores)
  }
  
  
  ExtractSegm <- as.data.frame(Reduce(rbind,seg.test))
  colnames(ExtractSegm) <- c("cell", "startseg","endseg","mean","sdVect")
  
  return(ExtractSegm)
}


getExtractSegm_old <- function(MatrixSeg){
  
  cellSegm <- MatrixSeg[,1]
  breaks <- which(diff(cellSegm)!=0)
  startseg <- c(1,(1+breaks))
  endseg <- c(breaks,length(cellSegm))
  
  getCellSegm <- function(i){
    cellSegm <- MatrixSeg[,i]
    sdVect <- rep.int(0,length(startseg))
    mean <- cellSegm[startseg]
    cell <- rep.int(i,length(startseg))
    cbind(cell, startseg, endseg, mean,sdVect)
  }
  
  ExtractSegm <-  as.data.frame(do.call(rbind,lapply(1:ncol(MatrixSeg), getCellSegm)))
  
  ExtractSegm
}


getLabelCall <- function(posteriorProbability){
  CallResults <- max.col(posteriorProbability) - 1
  
  ProbResults <- apply(posteriorProbability,1,max)
  
  #CallResults[ProbResults < 0.60] <- 0
  
  Results <- as.data.frame(cbind(CallResults,ProbResults))
  Results
}


getCallingCN <- function(AnnotMatrix,ExtractSegm, truncBoundRight,truncBoundLeft,meanVect, organism = "human", par_cores = 20){

  #meanVect <- c(-truncBoundLeft*2,-truncBoundLeft*1.5,0,truncBoundRight*1.5,truncBoundRight*2)
  #meanVect <- c(-truncBoundLeft*4,-truncBoundLeft*2,0,truncBoundRight*2,truncBoundRight*4)
  sdVect <- c(0.01,0.01,0.01,0.01,0.01)
  
  prior <- c(0.05,0.1,0.7,0.1,0.05)
  
  ResultsEM <- EM(ExtractSegm$mean,meanVect,sdVect,truncBoundRight,truncBoundLeft, prior)
  
  posteriorProbability <- getPosteriorProbability(ExtractSegm$mean,ResultsEM$meanVect,ResultsEM$sdVect,ResultsEM$prior)
  out <- getLabelCall(posteriorProbability)
  
  df_call <- as.data.frame(cbind(AnnotMatrix[ExtractSegm[,2],1],AnnotMatrix[ExtractSegm[,2],2],AnnotMatrix[ExtractSegm[,3],3],ExtractSegm[,4],out))
  colnames(df_call) <- c("Chromosome","Start","End","Mean","CN","ProbCall")
  
  getMajorityCall <- function(ch){
    call <- c()
    df_call_ch <- df_call[df_call$Chromosome==ch,]
    for(POSini in unique(df_call_ch$Start)){
      call <- c(call,as.numeric(names(which.max(table(df_call_ch[df_call_ch$Start==POSini,]$CN)))))
    }
    call
  }
  
  if(organism == "human"){
    totChr <- 22
  }else{
    totChr <- 19
  }   
  
  library(parallel)
  
  
  if(Sys.info()["sysname"]=="Windows"){
    cl <- parallel::makeCluster(getOption("cl.cores", par_cores))
    call <- parallel::parLapply(cl, 1:totChr, getMajorityCall)
    parallel::stopCluster(cl)
  }else{
    call <- parallel::mclapply(1:totChr, getMajorityCall ,mc.cores = par_cores)
  }
  
  call <- unlist(call)
  
  CNV <- df_call[1:which(df_call$Chromosome==totChr & df_call$End == max(df_call[df_call$Chromosome==totChr,]$End))[1],c(1,2,3,4)]
  CNV$Mean <- call
  colnames(CNV) <- c("Chr","Pos","End","CN")
  #plotSegmentation(CNV)
  
  CNV
}


getCNcall <- function(MatrixSeg, count_mtx_annot, breaks, sample = "",subclone = "", par_cores = 20, CLONAL = FALSE, organism = "human"){
  
  if(CLONAL){
    truncBoundLeft <- 0.05
    truncBoundRight <- 0.05
    meanVect <- c(-truncBoundLeft*4,-truncBoundLeft*2,0,truncBoundRight*2,truncBoundRight*4)
  }else{
    truncBoundLeft <- 0.1
    truncBoundRight <- 0.1
    meanVect <- c(-truncBoundLeft*2,-truncBoundLeft*1.5,0,truncBoundRight*1.5,truncBoundRight*2)
  }
  
  ExtractSegm <- getExtractSegm(MatrixSeg, breaks, par_cores)
  
  CNV <- getCallingCN(count_mtx_annot[,c(1,2,3)], ExtractSegm, truncBoundRight, truncBoundLeft, meanVect, organism = organism, par_cores = par_cores)
  #plotSegmentation(CNV)
  #write.table(CNV, file = paste("./output/",sample,"_",subclone,"_CN.seg"), sep = "\t", quote = FALSE)
  CNV
}

