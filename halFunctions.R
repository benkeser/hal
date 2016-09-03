#================================================
# Highly Adaptive Lasso and associated functions
#================================================
SL.hal <- function(Y, X, newX, family=gaussian(), 
                   verbose=TRUE,
                   obsWeights=rep(1,length(Y)), 
                   sparseMat=TRUE,
                   nfolds = ifelse(length(Y)<=100, 20, 10), 
                   nlambda = 100, useMin = TRUE,...){
  require(Matrix)
  d <- ncol(X)
  if(!sparseMat){
    uniList <- alply(as.matrix(X),2,function(x){
      # myX <- matrix(x,ncol=length(unique(x)), nrow=length(x)) - 
      #   matrix(unique(x), ncol=length(unique(x)), nrow=length(x), byrow=TRUE)
      myX <- matrix(x,ncol=length(x), nrow=length(x)) -
        matrix(x, ncol=length(x), nrow=length(x), byrow=TRUE)
      myX <- ifelse(myX < 0, 0, 1)
      myX
    })
    
    if(d >=2){
      highDList <- alply(matrix(2:d),1,function(k){
        thisList <- alply(combn(d,k),2,function(x){
          Reduce("*",uniList[x])
        })
        Reduce("cbind",thisList)
      })
      initX <- cbind(Reduce("cbind",uniList), Reduce("cbind",highDList))
      dup <- duplicated(t(initX))
      designX <- initX[,!dup]
    }else{
      initX <- Reduce("cbind",uniList)
      dup <- duplicated(t(initX))
      designX <- initX[,!dup]
    }
    
    fitCV <- glmnet::cv.glmnet(x = designX, y = Y, weights = obsWeights, 
                               lambda.min.ratio=0.001,
                               lambda = NULL, type.measure = "deviance", nfolds = nfolds, 
                               family = family$family, alpha = 1, nlambda = nlambda)
    
    
    ## get predictions back
    mynewX <- matrix(newX[,1],ncol=length(X[,1]), nrow=length(newX[,1])) - 
      matrix(X[,1], ncol=length(X[,1]), nrow=length(newX[,1]), byrow=TRUE)
    mynewX <- ifelse(mynewX < 0, 0, 1)
    
    makeNewDesignX <- TRUE
    if(all(dim(X)==dim(newX)))
      makeNewDesignX <- !all(X==newX)
    
    if(makeNewDesignX){
      uniList <- alply(matrix(1:ncol(X)),1,function(x){
        myX <- matrix(newX[,x],ncol=length(X[,x]), nrow=length(newX[,x])) - 
          matrix(X[,x], ncol=length(X[,x]), nrow=length(newX[,x]), byrow=TRUE)
        myX <- ifelse(myX < 0, 0, 1)
        myX
      })
      
      if(d >=2){
        highDList <- alply(matrix(2:d),1,function(k){
          thisList <- alply(combn(d,k),2,function(x){
            Reduce("*",uniList[x])
          })
          Reduce("cbind",thisList)
        })
        
        initX <- cbind(Reduce("cbind",uniList), Reduce("cbind",highDList))
        designNewX <- initX[,!dup]
      }else{
        initX <- Reduce("cbind",uniList)
        designNewX <- initX[,!dup]
      }
    }else{
      designNewX <- designX
    }
    
    pred <- predict(fitCV$glmnet.fit, newx = designNewX, 
                    s = ifelse(useMin,fitCV$lambda.min, fitCV$lambda.1se), type = "response")
    fit <- list(object = fitCV, useMin = useMin, X=X, dup=dup, sparseMat=sparseMat)
  }else{
    SuperLearner:::.SL.require("plyr")
    SuperLearner:::.SL.require("data.table")
    SuperLearner:::.SL.require("stringr")
    SuperLearner:::.SL.require("bit")
    
    if(is.vector(X)) X <- matrix(X, ncol=1)
    if(is.vector(newX)) newX <- matrix(newX, ncol=1)
    n <- length(X[,1])
    d <- ncol(X)
    
    if(verbose) cat("Making sparse matrix \n")
    X.init <- makeSparseMat(X=X,newX=X,verbose=TRUE)
    
    ## find duplicated columns
    if(verbose) cat("Finding duplicate columns \n")
    # Number of columns will become the new number of observations in the data.table
    nIndCols <- ncol(X.init)
    # Pre-allocate a data.table with one column, each row will store a single column from X.init
    datDT <- data.table(ID = 1:nIndCols, bit_to_int_to_str = rep.int("0", nIndCols))
    # Each column in X.init will be represented by a unique vector of integers.
    # Each indicator column in X.init will be converted to a row of integers or a string of cat'ed integers in data.table
    # The number of integers needed to represent a single column is determined automatically by package "bit" and it depends on nrow(X.init)
    nbits <- nrow(X.init) # number of bits (0/1) used by each column in X.init
    bitvals <- bit(length = nbits) # initial allocation (all 0/FALSE)
    nints_used <- length(unclass(bitvals)) # number of integers needed to represent each column
    # For loop over columns of X.init
    ID_withNA <- NULL # track which results gave NA in one of the integers
    for (i in 1:nIndCols) {
      bitvals <- bit(length = nbits) # initial allocation (all 0/FALSE)
      Fidx_base0 <- (X.init@p[i]) : (X.init@p[i + 1]-1) # zero-base indices of indices of non-zero rows for column i=1
      nonzero_rows <- X.init@i[Fidx_base0 + 1] + 1 # actual row numbers of non-zero elements in column i=1
      # print(i); print(nonzero_rows)
      # X.init@i[i:X.init@p[i]]+1 # row numbers of non-zero elements in first col
      bitvals[nonzero_rows] <- TRUE
      # str(bitwhich(bitvals))
      intval <- unclass(bitvals) # integer representation of the bit vector
      # stringval <- str_c(intval, collapse = "")
      if (any(is.na(intval))) ID_withNA <- c(ID_withNA, i)
      set(datDT, i, 2L, value = str_c(str_replace_na(intval), collapse = ""))
    }
    # create a hash-key on the string representation of the column, 
    # sorts it by bit_to_int_to_str using radix sort:
    setkey(datDT, bit_to_int_to_str)
    # add logical column indicating duplicates, 
    # following the first non-duplicate element
    datDT[, duplicates := duplicated(datDT)]
    # just get the column IDs and duplicate indicators:
    datDT[, .(ID, duplicates)]
    
    dupInds <- datDT[,ID][which(datDT[,duplicates])]
    uniqDup <- unique(datDT[duplicates==TRUE,bit_to_int_to_str])
    
    colDups <- alply(uniqDup, 1, function(l){
      datDT[,ID][which(datDT[,bit_to_int_to_str] == l)]
    })
    
    if(verbose) cat("Fitting lasso \n")
    if(length(dupInds)>0){
      notDupInds <- (1:ncol(X.init))[-unlist(colDups, use.names = FALSE)]
      keepDupInds <- unlist(lapply(colDups, function(x){ x[[1]] }), use.names=FALSE)
      
      fitCV <- glmnet::cv.glmnet(x = X.init[,c(keepDupInds,notDupInds)], y = Y, weights = obsWeights, 
                                 lambda = NULL, lambda.min.ratio=0.001, type.measure = "deviance", nfolds = nfolds, 
                                 family = family$family, alpha = 1, nlambda = nlambda)
    }else{
      fitCV <- glmnet::cv.glmnet(x = X.init, y = Y, weights = obsWeights, 
                                 lambda = NULL, lambda.min.ratio=0.001, type.measure = "deviance", nfolds = nfolds, 
                                 family = family$family, alpha = 1, nlambda = nlambda)
    }
    
    fit <- list(object = fitCV, useMin = useMin, X=X, dupInds = dupInds, colDups=colDups, sparseMat=sparseMat)
    class(fit) <- "SL.hal"
    
    if(identical(X,newX)){
      if(length(dupInds) > 0){
        pred <- predict(fitCV, newx = X.init[,c(keepDupInds,notDupInds)], s = ifelse(useMin, fitCV$lambda.min, fitCV$lambda.1se), 
                        type = "response")
      }else{
        pred <- predict(fitCV, newx = X.init, s = ifelse(useMin, fitCV$lambda.min, fitCV$lambda.1se), 
                        type = "response")
      }
    }else{
      pred <- predict(fit, newdata=newX, bigDesign=FALSE, chunks=10000)          
    }
  }
  
  out <- list(pred = pred, fit = fit)
  cat("Done with SL.hal")
  return(out)
}



predict.SL.hal <- function (object, newdata, bigDesign=FALSE, verbose=TRUE, 
                            chunks=1000,
                            s = ifelse(object$useMin, object$object$lambda.min, object$object$lambda.1se),...)
{
  if(!object$sparseMat){
    d <- ncol(object$X)
    # if you want to get predictions all at once (smaller newdata)
    if(bigDesign){
      uniList <- alply(matrix(1:ncol(object$X)),1,function(x){
        myX <- matrix(newdata[,x],ncol=length(object$X[,x]), nrow=length(newdata[,x])) - 
          matrix(object$X[,x], ncol=length(object$X[,x]), nrow=length(newdata[,x]), byrow=TRUE)
        myX <- ifelse(myX < 0, 0, 1)
        myX
      })
      
      if(d >=2){
        highDList <- alply(matrix(2:d),1,function(k){
          thisList <- alply(combn(d,k),2,function(x){
            Reduce("*",uniList[x])
          })
          Reduce("cbind",thisList)
        })
        
        initX <- cbind(Reduce("cbind",uniList), Reduce("cbind",highDList))
        designNewX <- initX[,!object$dup]
      }else{
        initX <- Reduce("cbind",uniList)
        designNewX <- initX[,!object$dup]
      }
      # get predictions
      pred <- predict(object$object$glmnet.fit, newx = designNewX, 
                      s = s, 
                      type = "response")
      
    }else{
      # get row by row predictions, so you never have to store a big design matrix
      # for newdata
      pred <- apply(as.matrix(newdata),1,function(i){
        uniList <- alply(matrix(1:ncol(object$X)),1,function(x){
          myX <- matrix(i[x],ncol=length(object$X[,x]), nrow=length(i[x])) - 
            matrix(object$X[,x], ncol=length(object$X[,x]), nrow=length(i[x]), byrow=TRUE)
          myX <- ifelse(myX < 0, 0, 1)
          myX
        })
        
        if(d >=2){
          highDList <- alply(matrix(2:d),1,function(k){
            thisList <- alply(combn(d,k),2,function(x){
              Reduce("*",uniList[x])
            })
            Reduce("cbind",thisList)
          })
          
          initX <- cbind(Reduce("cbind",uniList), Reduce("cbind",highDList))
          designNewX <- initX[!object$dup]
        }else{
          initX <- Reduce("cbind",uniList)
          designNewX <- initX[!object$dup]
        }
        # get predictions
        thispred <- predict(object$object$glmnet.fit, newx = matrix(designNewX,nrow=1), s=s,
                            type = "response")
        thispred
      })
    }
  }else{
    # all predictions at once
    if(bigDesign){
      pred <- doPred(object=object,newdata=newdata,verbose=verbose)
    }else{
      nNew <- length(newdata[,1])
      nChunks <- floor(nNew/chunks) + ifelse(nNew%%chunks==0, 0, 1)
      pred <- rep(NA, length(nNew))
      for(i in 1:nChunks){
        minC <- (i-1)*chunks + 1
        maxC <- ifelse(i==nChunks, nNew,i*chunks)
        pred[minC:maxC] <- doPred(object=object,s=s,newdata=newdata[minC:maxC,],verbose=verbose)
      }
    }
  }
  return(as.numeric(pred))
}

doPred <- function(object,newdata,verbose=TRUE,s){
  if(is.vector(newdata)) newdata <- matrix(newdata)
  
  if(verbose) cat("Making initial sparse matrix for predictions \n")
  tmp <- makeSparseMat(X=object$X, newX=newdata, verbose=verbose)
  
  if(length(object$dupInds) > 0){
    if(verbose) cat("Correcting duplicate columns in sparse matrix \n")
    # get vector of duplicate columns
    dupVec <- unlist(object$colDups,use.names=FALSE)
    # number of each duplicate
    nperDup <- unlist(lapply(object$colDups,length),use.names=FALSE)
    # number of duplicates to roll through
    K <- length(nperDup)
    # start and ending index
    startInd <- c(0, cumsum(nperDup)[1:(K-1)])
    endInd <- cumsum(nperDup)
    # not duplicate colums
    notdupVec <- (1:ncol(tmp))[-dupVec]
    # put all the duplicated guys first
    tmp <- tmp[,c(dupVec,notdupVec)]
    
    uniqRowsList <- list()
    myp <- c(0,rep(NA, K))
    # look at the i associatiated with 
    for(k in 1:K){
      # this condition ensures that not all the values of a given set of duplciates
      # are equal to zero.
      if(tmp@p[startInd[k]+1] != tmp@p[endInd[k]+1]){
        Fidx_base0 <- (tmp@p[startInd[k] + 1]) : (tmp@p[endInd[k] + 1] - 1)
        nonzero_rows <- tmp@i[Fidx_base0 + 1] + 1 # actual row numbers of non-zero elements in column i=1
        # unique nonzero_rows
        theseRows <- sort(unique(nonzero_rows))
        uniqRowsList[[k]] <- theseRows
        # a new p vector for duplicated columns
        myp[k+1] <- myp[k] + length(theseRows)
      }else{ # the whole block for this set of duplicates is 0
        uniqRowsList[[k]] <- NULL
        myp[k+1] <- myp[k]
      }
    }
    
    # look at this sparse matrix
    myi <- unlist(uniqRowsList)
    # check if it came out right
    # grbg1 <- sparseMatrix(i=myi, p=myp, x=1)
    
    # fix up p with nondup columns
    ## for this example every non-duplicated column in the new design
    ## matrix is 0, which is causing this to break. I think. 
    if(tmp@p[endInd[K] + 1] != tmp@p[length(tmp@p)]){
      fulli <- c(myi, tmp@i[(tmp@p[endInd[K] + 1] + 1) : tmp@p[length(tmp@p)]] + 1)
      fullp <- c(myp, 
                 tmp@p[((endInd[K]+1) + 1) : length(tmp@p)] - 
                   tmp@p[(endInd[K]+1)] + myp[K+1])
    }else{
      fulli <- myi
      fullp <- myp
    }
    # 
    tmp <- sparseMatrix(i=fulli, p=fullp, x=1, 
                        dims=c(length(newdata[,1]),
                               length(notdupVec) + length(object$colDup)))
  }
  pred <- predict(object$object$glmnet.fit, newx=tmp, 
                  s = s)
  pred
}