## Cell type prediction from minfi adapted by Eilis Hannon
## functions for estimate cell counts edited to take matrices rather than mset
## Error metric by Dorothea Seiler Vellame

library(genefilter)
library(quadprog)


projectCellTypeWithError = function(YIN, modelType = c("mouseHumanHybrid", "mouse", "ownModel"), ownModelData = NA){
  
  sampleDup = 0
  if(modelType != "ownModel"){
    load(paste("models/", modelType, "Model.Rdata", sep = ""))
  }else{
    model = ownModelData
  }
  
  ## subset the CpGs in YIN 
  YIN = YIN[rownames(model$coefEsts),]
  coefCellTypeIN = model$coefEsts
  
  ## if there's only 1 sample, pretend there are 2 so it doesn't break
  if (ncol(as.matrix(YIN)) == 1){ 
    sampleDup = 1
    YIN = cbind(YIN, YIN)
  }
  
  projectCellType <- function(Y, coefCellType, contrastCellType=NULL, nonnegative=TRUE, lessThanOne=FALSE){ 
    if(is.null(contrastCellType))
      Xmat <- coefCellType
    else
      Xmat <- tcrossprod(coefCellType, contrastCellType) 
    
    nCol <- dim(Xmat)[2]
    if(nCol == 2) {
      obs <- which(apply(Y, 1, function(x){return(sum(is.na(x)) == 0)}))
      Dmat <- crossprod(Xmat[obs,])
      mixCoef <- t(apply(Y, 2, function(x) { solve(Dmat, crossprod(Xmat[obs,], x[obs])) }))
      colnames(mixCoef) <- colnames(Xmat)
      return(mixCoef)
    } else {
      nSubj <- dim(Y)[2]
      
      mixCoef <- matrix(0, nSubj, nCol)
      rownames(mixCoef) <- colnames(Y)
      colnames(mixCoef) <- colnames(Xmat)
      
      if(nonnegative){
        if(lessThanOne) {
          Amat <- cbind(rep(-1, nCol), diag(nCol))
          b0vec <- c(-1, rep(0, nCol))
        } else {
          Amat <- diag(nCol)
          b0vec <- rep(0, nCol)
        }
        for(i in 1:nSubj) {
          obs <- which(!is.na(Y[,i])) 
          Dmat <- crossprod(Xmat[obs,])
          mixCoef[i,] <- solve.QP(Dmat, crossprod(Xmat[obs,], Y[obs,i]), Amat, b0vec)$sol
        }
      } else {
        for(i in 1:nSubj) {
          obs <- which(!is.na(Y[,i])) 
          Dmat <- crossprod(Xmat[obs,])
          mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
        }
      }
      return(mixCoef)
    }
  }
  
  ## error function
  getErrorPerSample = function(applyIndex,
                               predictedIN = mixCoef,
                               coefDataIN = coefCellTypeIN,
                               betasBulkIN = YIN){
    
    trueBulk = matrix(ncol = 1, nrow = nrow(coefDataIN), data = 0)
    
    RMSE = function(m, o){
      sqrt(mean((m - o)^2))
    }
    
    for (i in 1:ncol(coefDataIN)){
      
      trueBulk[,1] = trueBulk[,1] + coefDataIN[,i]*predictedIN[applyIndex,i]
    }
    
    betasBulkIN = t(apply(betasBulkIN, 1, function(x){x[is.na(x)] = 0; return(x)}))
    
    error = RMSE(trueBulk, betasBulkIN[,applyIndex])
    return(error)
  }
  
  mixCoef = projectCellType(YIN, coefCellTypeIN)
  error = sapply(1:nrow(mixCoef), getErrorPerSample)
  nCGmissing = apply(YIN, 2, function(x){sum(is.na(x))})
  mixCoef = cbind(mixCoef, error, nCGmissing)
  
  if (sampleDup == 1){
    mixCoef = mixCoef[1,]
  }
  return(mixCoef)
}