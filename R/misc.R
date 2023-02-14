## non-exported functions from MTAGEI package
recoverLD <- function(sumstats, snp.list, bymean = TRUE){
  ## recover R.C from input of V, so MTAGEI doesn't need the extra input of R.C
  ## start ##
  D <- length(sumstats)
  K <- length(sumstats[[1]])

  R.C <- matrix(0, length(snp.list), length(snp.list))
  diag(R.C) <- 1
  colnames(R.C) <- rownames(R.C) <- snp.list
  row_col.id <- which(lower.tri(R.C), arr.ind = TRUE)

  R.bydrugs <- list()
  for(d in 1:D) {
    R.bytraits <- list()
    for(k in 1:K){
      R.bytraits[[k]] <- stats::cov2cor(sumstats[[d]][[k]]$V)
    }
    R.bydrugs[[d]] <- R.bytraits
  }

  for (iter in 1:nrow(row_col.id)) {
    row_col.name <- c(rownames(R.C)[row_col.id[iter, 1]],
                      colnames(R.C)[row_col.id[iter, 2]])
    # stop <- FALSE
    R.tmp <- NULL
    for (d in 1:D) {
      for (k in 1:K) {
        R.bytraits <- R.bydrugs[[d]][[k]]
        row.id <- which(rownames(R.bytraits) == row_col.name[1])
        col.id <- which(colnames(R.bytraits) == row_col.name[2])
        if (length(row.id) != 0 & length(col.id) != 0) {
          if (!is.nan(R.bytraits[row.id, col.id])) {
            R.tmp <- c(R.tmp, R.bytraits[row.id, col.id])
          }
        }
      }
    }
    if(is.null(R.tmp)) {
      R.C[row_col.id[iter, 1], row_col.id[iter, 2]] <- NA
    }else{
      if(bymean) {
        R.C[row_col.id[iter, 1], row_col.id[iter, 2]] <- mean(R.tmp, na.rm = T)
      }else{
        R.C[row_col.id[iter, 1], row_col.id[iter, 2]] <- max(R.tmp, na.rm = T)

      }
    }
  }
  R.C <- Matrix::forceSymmetric(R.C, uplo = "L")
  R.C <- as.matrix(R.C)
  R.C[is.na(R.C)] <- 0 #this would happen if there is no/partial overlap and different traits have different list of nonpolymorphic SNPs
  return(R.C)
}
nonpolymorphic.fn <- function(mat, obs.stat, GEI = FALSE, nonpoly.list = NULL) {
  ## remove nonpolymorphic SNPs in calculating GEI effects
  D <- length(obs.stat)
  K <- length(obs.stat[[1]])
  trait.ind <- expand.grid(row.ind = 1:K, column.ind = 1:K)

  mat.bydrug <- list()
  for(d in 1:D) {
    mat1 <- list()
    for(iter in 1:nrow(trait.ind)) {
      k1 <- trait.ind[iter, 1]
      k2 <- trait.ind[iter, 2]
      if(!GEI) {
        ind1 <- which(obs.stat[[d]][[k1]]$U == 0)
        ind2 <- which(obs.stat[[d]][[k2]]$U == 0)
      }else{
        ind1 <- nonpoly.list[[k1]]
        ind2 <- nonpoly.list[[k2]]
      }

      if(length(ind1) != 0 & length(ind2) != 0) {
        mat1[[iter]] <- mat[-ind1, -ind2, drop = FALSE]
      } else if(length(ind1) != 0 & length(ind2) == 0) {
        mat1[[iter]] <- mat[-ind1, , drop = FALSE]
      } else if(length(ind1) == 0 & length(ind2) != 0) {
        mat1[[iter]] <- mat[, -ind2, drop = FALSE]
      } else {
        mat1[[iter]] <- mat
      }
    }
    mat.bydrug[[d]] <- mat1
  }

  return(mat.bydrug)
}

remove.nonpoly <- function(x, nonpoly.k1, nonpoly.k2) {
  if(inherits(x, "numeric")) { #if (class(x) == "numeric") { a bad form to check that the class was equal to a particular value
    if(length(nonpoly.k1) != 0) {
      ind1 <- which(names(x) %in% nonpoly.k1)
      if(length(ind1) != 0) {
        x <- x[-ind1]
      }
    }
  }else if(inherits(x, "matrix")) {#if(class(x) == "matrix") {
    if(length(nonpoly.k1)!=0 & length(nonpoly.k2) != 0) {
      ind1 <- which(rownames(x) %in% nonpoly.k1)
      ind2 <- which(colnames(x) %in% nonpoly.k2)
      if(length(ind1) != 0) {
        x <- x[-ind1, , drop = FALSE]
      }
      if(length(ind2) != 0) {
        x <- x[, -ind2, drop = FALSE]
      }

    }else if(length(nonpoly.k1)==0 & length(nonpoly.k2) != 0) {
      ind2 <- which(colnames(x) %in% nonpoly.k2)
      if(length(ind2) != 0) {
        x <- x[, -ind2, drop = FALSE]
      }
    }else if(length(nonpoly.k1) !=0 & length(nonpoly.k2) == 0) {
      ind1 <- which(rownames(x) %in% nonpoly.k1)
      if(length(ind1) != 0) {
        x <- x[-ind1, , drop = FALSE]
      }
    }
  }
  return(x)
}
count.overlapping <- function(overlapping.samples){
  K <- length(overlapping.samples)
  trait.grid <- gtools::combinations(K, 2, 1:K, repeats.allowed = TRUE)
  trait.count <- matrix(0, K, K)
  for(row.id in 1:nrow(trait.grid)){
    k1 <- trait.grid[row.id, 1]
    k2 <- trait.grid[row.id, 2]
    trait.count[k1, k2] <- length(intersect(overlapping.samples[[k1]], overlapping.samples[[k2]]))
  }
  trait.count <- trait.count + t(trait.count) - diag(diag(trait.count), nrow = K, ncol = K)
  return(trait.count)
}

SPA.convertV.bygrp <- function(SNP.name, d, beta.sumstats, beta.sumstats.SPA, K){
  V.origin <- matrix(sapply(beta.sumstats[[d]]$beta.cov, function(x)
    x[rownames(x) == SNP.name, rownames(x) == SNP.name]), K, K, byrow = FALSE)
  R.origin <- cov2cor(V.origin)
  V.SPA.sqrt.inv <- diag(sqrt(1/sapply(beta.sumstats.SPA[[d]]$V.SPA, function(x) x[names(x) == SNP.name])))
  cov.SPA <- V.SPA.sqrt.inv %*% R.origin %*% V.SPA.sqrt.inv
  trait.grid <- expand.grid(trait1 = 1:K, trait2 = 1:K)

  ret <- list()
  for(i in 1:nrow(trait.grid)) {
    temp <- cov.SPA[trait.grid[i, 1], trait.grid[i, 2], drop = FALSE]
    rownames(temp) <- colnames(temp) <- SNP.name
    ret[[i]] <- temp
  }
  names(ret) <- apply(trait.grid, 1, function(x) paste0("Traits", x[1], "_", x[2]))
  return(ret)
}

Beta.Weights<-function(MAF, weights.beta = c(1, 25), Cutoff=1, Is.MAF=TRUE){
  n<-length(MAF)
  weights<-rep(0,n)
  Sign<-rep(1,n)
  if(!is.null(weights.beta)){
    IDX1<-which(MAF > 0.5)
    if(length(IDX1) > 0){
      Sign[IDX1]<--1
      MAF[IDX1]<-1-MAF[IDX1]
    }

    IDX_0<-union(which(MAF == 0), which(MAF > Cutoff))
    if(length(IDX_0) == n){
      #stop("No polymorphic SNPs")
      weights<-rep(0,n)
    } else if( length(IDX_0) == 0){
      weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
    } else {
      weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
    }

    weights = weights * Sign
  }else{
    weights <- rep(1, n)
  }
  return(weights)

}

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  rho^exponent
}
combMTAR <- function(x, rho){
  y <- sign(x) * abs(x)^rho
  return(y)
}
