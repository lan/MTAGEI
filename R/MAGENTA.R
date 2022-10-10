ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  rho^exponent
}
combMTAR <- function(x, rho){
  y <- sign(x) * abs(x)^rho
  return(y)
}

MAGENTA <- function(beta.sumstats, MAF, input = "data",
                    test = c("joint", "GEI"),
                    MAF.thres,
                    KA, zeta.ret,
                    drugStruct = "AR1",
                    # diagTrans = FALSE,
                    R.C = NULL, cct = TRUE,
                    diffweight = FALSE,
                    threshold = 0.05,
                    weight.commonSNP = NULL,
                    pval.thres = 1,
                    rho.SNP = c(0, 0.5, 1),
                    rho.trait = c(0.5, 1.5, 1),
                    rho.drug = c(0, 0.5, 1),
                    weight.SNP = c(1, 25)){
  # print(paste0("Start analysis at ", Sys.time()))
  D <- length(sumstats)
  K <- length(sumstats[[1]])
  m <- length(MAF)
  snp.list <- names(MAF)
  drug_trait <- expand.grid(trait = 1:K, drug = 1:D)

  zeta.GE <- matrix(0, D*K, D*K)
  for(d in 1:D) {
    ind <- 1:K + (d - 1) * K
    zeta.GE[ind, ind] <- zeta.list[[d]]
  }
  zeta.GE1 <- diag(D)

  # recalculate V according to R.C
  for(d in 1:D) {
    sumstats.bydrug <- list()
    for(k in 1:K) {
      if(is.matrix(V)) {
        V.sqrt.tmp <- diag(nrow(sumstats[[d]][[k]]$V))
        diag(V.sqrt.tmp) <- sqrt(diag(sumstats[[d]][[k]]$V))
      }else{
        V.sqrt.tmp <- diag(m)
        diag(V.sqrt.tmp) <- sqrt(sumstats[[d]][[k]]$V)
      }

      V.tmp <- V.sqrt.tmp %*% R.C %*% V.sqrt.tmp
      rownames(V.tmp) <- colnames(V.tmp) <- colnames(sumstats[[d]][[k]]$V)
      sumstats.bydrug[[k]] <- list(U = sumstats[[d]][[k]]$U,
                                   V = V.tmp)
    }
    sumstats[[d]] <- sumstats.bydrug
  }

  U.comb <- V.comb <- list()
  for(k in 1:K) {
    U.temp <- sumstats[[1]][[k]]$U
    V.temp <- sumstats[[1]][[k]]$V
    for(d in 2:D) {
      U.temp <- U.temp + sumstats[[d]][[k]]$U
      V.temp <- V.temp + sumstats[[d]][[k]]$V
    }
    U.comb[[k]] <- U.temp
    V.comb[[k]] <- V.temp
  }

  sumstats.bytrait <- list()
  zeta.list1 <- list()
  for(k in 1:K) {
    sumstats.bytrait.tmp <- list()
    for(d in 1:D) {
      zeta.list1[[d]] <- matrix(1)
      sumstats.bytrait.tmp[[d]] <- list(list(U = sumstats[[d]][[k]]$U, V = sumstats[[d]][[k]]$V))
    }
    sumstats.bytrait[[k]] <- sumstats.bytrait.tmp
  }

  ## main effect test start ##
  combMTAR.p <- try(MAGENTA.main(U = U.comb, V = V.comb, MAF = MAF, R.C = R.C,
                                 snp.list = names(MAF),
                                 # diagTrans = diagTrans,
                                 cct = cct,
                                 KA = KA, zeta = zeta, diffweight = diffweight)$p, silent = T)
  singletrait.p <- list()
  for(k in 1:K){
    singletrait.p[[k]] <- try(MAGENTA.main(U = list(U.comb[[k]]), V = list(V.comb[[k]]), MAF = MAF, R.C = R.C,
                                           snp.list = names(MAF),
                                           # diagTrans = diagTrans,
                                           KA = matrix(1), rho.trait = 1,
                                           zeta = matrix(1), cct = cct,
                                           diffweight = diffweight)$p,silent = T)
  }

  if(inherits(combMTAR.p, "try-error") | any(sapply(singletrait.p, function(x)
    inherits(x, "try-error")))) {
    main.p <- c(NA, rep(NA, K))
  }else{
    main.p <- c(combMTAR.p, unlist(singletrait.p))

  }
  names(main.p) <- c("MTARC.M", paste("Trait", 1:K, "M", sep = "."))


  # joint effect
  combMTAR.p <- try(MAGENTA.joint(sumstats = sumstats, MAF = MAF, KA = KA,
                                  # diagTrans = diagTrans,
                                  R.C = R.C, zeta = zeta.list,
                                  drugStruct = drugStruct,
                                  rho.drug = c(0, 0.5, 1), cct = cct,
                                  diffweight = diffweight)$p,silent = T)

  singletrait.p <- list()
  singletrait.time <- rep(NA, K)
  for(k in 1:K){
    singletrait.p[[k]] <- try(MAGENTA.joint(sumstats = sumstats.bytrait[[k]],
                                            MAF = MAF, KA = matrix(1),
                                            # diagTrans = diagTrans,
                                            R.C = R.C, zeta = zeta.list1,
                                            drugStruct = drugStruct,
                                            rho.drug = c(0, 0.5, 1), cct = cct,
                                            rho.trait = 1,
                                            diffweight = diffweight)$p, silent = T)

  }

  if(inherits(combMTAR.p, "try-error") |
     any(sapply(singletrait.p, function(x) inherits(x, "try-error")))) {
    joint.p <- c(NA, rep(NA, K))
  }else{
    joint.p <- c(combMTAR.p, unlist(singletrait.p))
  }
  names(joint.p) <- c("MTARC.J", paste("Trait", 1:K, "J", sep = "."))

  ## interaction effect test start ##
  combMTAR.p <- try(MAGENTA.GEI(sumstats = sumstats, MAF = MAF, KA = KA,
                                # diagTrans = diagTrans,
                                R.C = R.C, zeta = zeta.GE, drugStruct = drugStruct,
                                rho.drug = c(0, 0.5, 1), cct = cct,
                                diffweight = diffweight)$p,silent = T)

  singletrait.p <- list()
  for(k in 1:K){
    singletrait.p[[k]] <- try(MAGENTA.GEI(sumstats = sumstats.bytrait[[k]], MAF = MAF,
                                          KA = matrix(1),
                                          # diagTrans = diagTrans,
                                          R.C = R.C,
                                          zeta = zeta.GE1, drugStruct = drugStruct,
                                          rho.drug = c(0, 0.5, 1), cct = cct,
                                          rho.trait = 1, diffweight = diffweight)$p, silent = T)

  }
  if(inherits(combMTAR.p, "try-error") |
     any(sapply(singletrait.p, function(x) inherits(x, "try-error")))) {
    interaction.p1 <- c(NA, rep(NA, K))
  }else{
    interaction.p1 <- c(combMTAR.p, unlist(singletrait.p))

  }
  names(interaction.p1) <- c("MTARC.IB", paste("Trait", 1:K, "IB", sep = "."))

  p.main.wo <- main.p
  p.joint.wo <- joint.p
  p.GEI.wo <- interaction.p1

  ## variantP ##
  MAC10 <- MAF[MAF > MAF.thres]
  highLD <- try(caret::findCorrelation(R.C[colnames(R.C) %in% names(MAC10),
                                           colnames(R.C) %in% names(MAC10)],
                                       cutoff = 0.98),silent = T)
  if(length(highLD)!=0 & !inherits(highLD, "try-error")) {
    selSNP <- names(MAC10)[-highLD]
  }else{
    selSNP <- names(MAC10)
  }

  ## way 3, wSPA, U = U, V = (U/Z)^2
  variant.p <- SPA_MAGENTA_diffU(simdata = simdata, genotype = genotype,
                                 sumstats = sumstats,
                                 zeta = zeta, zeta.list = zeta.list,
                                 zeta.GE = zeta.GE,
                                 selSNP = selSNP,
                                 MAF = MAF, R.C = R.C, KA = KA,
                                 MAC.thres = MAC.thres,way = 3)

  main.varp <- cbind(as.numeric(str_split(variant.p$main[1], ":")[[1]]),
                     matrix(as.numeric(str_split(variant.p$main[2], ":")[[1]]), ncol = 3))
  joint.varp <- cbind(as.numeric(str_split(variant.p$joint[1], ":")[[1]]),
                      matrix(as.numeric(str_split(variant.p$joint[2], ":")[[1]]), ncol = 3))
  GEI.varp <- cbind(as.numeric(str_split(variant.p$GEI[1], ":")[[1]]),
                    matrix(as.numeric(str_split(variant.p$GEI[2], ":")[[1]]), ncol = 3))
  rownames(main.varp) <- rownames(joint.varp) <- rownames(GEI.varp) <- selSNP

  p.main.w5 <- ACAT(c(p.main.wo, main.varp))
  p.joint.w5 <- ACAT(c(p.joint.wo, joint.varp))
  p.GEI.w5 <- ACAT(c(p.GEI.wo, GEI.varp))

  diffweight.p <- c(p.main.w5 = p.main.w5,
                    p.joint.w5 = p.joint.w5,
                    p.GEI.w5 = p.GEI.w5)

  tmp.ret <- data.frame(gene.size = m, MAC10 = length(which(MAF > MAF.thres)))

  tmp.ret <- cbind(tmp.ret, t(diffweight.p),
                   t(p.main.wo), t(p.joint.wo),
                   t(p.GEI.wo))
}

MAGENTA.joint <- function(sumstats, MAF, KA = NULL, drugStruct = NULL,
                          # diagTrans = FALSE,
                          R.C = NULL, cct = TRUE,
                          diffweight = FALSE,
                          threshold = 0.05,
                          weight.commonSNP = NULL,
                          zeta = NULL, pval.thres = 1,
                          rho.SNP = c(0, 0.5, 1),
                          rho.trait = c(0.5, 5, 1),
                          rho.drug = c(0, 0.5, 1),
                          weight.SNP = c(1, 25))
{
  # print(paste0("Start analysis at ", Sys.time()))
  D <- length(sumstats)
  K <- length(sumstats[[1]])
  m <- length(MAF)
  snp.list <- names(MAF)
  drug_trait <- expand.grid(trait = 1:K, drug = 1:D)

  if (is.null(drugStruct)){
    message("Assume the effects among environmental groups are independent")
    lambda <- expand.grid(lambdaA = rho.trait, lambdaC = rho.SNP, lambdaD = 0)
  }else if(drugStruct == "AR1") {
    message("Assume the first-order autoregressive structure for the effects among environmental groups")
    lambda <- expand.grid(lambdaA = rho.trait, lambdaC = rho.SNP, lambdaD = rho.drug)
  }

  if (is.null(KA)) {
    message("Without the input of genetic correlation information, an exchangeable correlation structure among traits is assumed in MTAR.")
    KA <- diag(K)
  }
  if (is.null(zeta)) {
    message("Without the input of zeta, MAGENTA assumes there are no overlap samples in the dataset.")
    zeta <- list()
    for(d in 1:D) {
      zeta[[d]] <- diag(K)
    }
  }
  message("Conducting MAGENTA analysis ...")

  if(diffweight) {
    message("Use different weights for rare and common variants:")
    if(is.null(weight.commonSNP)){
      message(paste0("No input of weight for common variants, use the default weight of dbeta(MAF,0.5,0.5) for common variants and dbeta(MAF,", paste0(weight.SNP, collapse = ","),  ") for rare variants"))
      weight.commonSNP <- c(0.5, 0.5)
    }else{
      message(paste0("Use the weight of dbeta(MAF,", paste0(weight.SNP, collapse = ","),  ") for rare variants and the weight of dbeta(MAF,", paste0(weight.commonSNP, collapse = ","), ") for common variants"))
    }
    common.SNP <- which(MAF >= threshold)
  }
  #
  #   # recover R.C from input of V, so MAGENTA doesn't need the extra input of R.C
  #   if(is.null(R.C)) {
  #     R.C <- recoverLD(sumstats = sumstats, snp.list = snp.list)
  #     message("The LD matrix is not provided and recovered from the summary statistics V.")
  #   }

  obs.stat <- list()
  U.complete <- list()
  for (d in 1:D) {
    U.temp <- list()
    obs.stat.temp <- list()
    for (k in 1:K) {
      U.bytrait <- numeric(length(snp.list))
      V.bytrait <- matrix(0, length(snp.list), length(snp.list))
      names(U.bytrait) <- snp.list
      colnames(V.bytrait) <- rownames(V.bytrait) <- snp.list
      order1 <- order(match(names(sumstats[[d]][[k]]$U)[names(sumstats[[d]][[k]]$U) %in% snp.list], snp.list))
      order2 <- which(snp.list %in% names(sumstats[[d]][[k]]$U))
      U.bytrait[order2] <- sumstats[[d]][[k]]$U[order1]
      V.bytrait[order2, order2] <- sumstats[[d]][[k]]$V[order1, order1]
      obs.stat.temp[[k]] <- list(U = U.bytrait, V = V.bytrait)
      U.temp[[k]] <- U.bytrait
    }
    obs.stat[[d]] <- obs.stat.temp
    U.complete[[d]] <- U.temp
  }


  obs.stat1 <- obs.stat

  if (!is.null(zeta)) {
    V.sqrt <- list() #contain the square root of all D by K covariance matrix V.
    for(d in 1:D) {
      V.sqrt.temp <- list()
      for(k in 1:K) {
        V.sqrt.temp[[k]] <- diag(length(obs.stat1[[d]][[k]]$U))
        diag(V.sqrt.temp[[k]]) <- sqrt(diag(obs.stat1[[d]][[k]]$V))
      }
      V.sqrt[[d]] <- V.sqrt.temp
    }

    trait.ind <- expand.grid(row.ind = 1:K, column.ind = 1:K)
    U.cov.bydrug <- list()
    V.diag.bydrug <- list()
    for(d in 1:D) {
      Ucov <- list()
      for (iter in 1:nrow(trait.ind)) {
        k1 <- trait.ind[iter, 1]
        k2 <- trait.ind[iter, 2]

        Ucov[[iter]] <- diag(V.sqrt[[d]][[k1]])  %*%  t(diag(V.sqrt[[d]][[k2]])) * R.C

        # if(!diagTrans) {
        #   if(k1 == k2) { #diagonal matrices
        #     Ucov[[iter]] <- obs.stat1[[d]][[k1]]$V
        #   }
        # }
      }
      U.cov <- NULL
      mlist <- list()
      for (col_id in 1:K) {
        diag.ind <- which(trait.ind$column.ind == col_id & trait.ind$row.ind == col_id)
        U.cov.col <- NULL
        index <- which(trait.ind$column.ind == col_id)
        row_id <- trait.ind[index, 1]
        for (iter in 1:length(index)) {
          U.cov.col <- rbind(U.cov.col, as.matrix(zeta[[d]][row_id[iter], col_id] * Ucov[[index[iter]]]))
        }
        U.cov <- cbind(U.cov, U.cov.col)
        mlist[[col_id]] <- Ucov[[diag.ind]]
      }
      U.cov <- as.matrix(U.cov)
      U.cov.bydrug[[d]] <- U.cov
      V.diag.bydrug[[d]] <- as.matrix(Matrix::bdiag(mlist))
    }
    U.cov.all <- as.matrix(Matrix::bdiag(U.cov.bydrug))

    U.inv <- try(MASS::ginv(U.cov.all), silent = T)
    # print(paste0("Done inverse of U at ", Sys.time()))
    if (inherits(U.inv, "try-error")) {
      warning("The covariance matrix of U is exactly singular, MASS:ginv() function doesn't work here.")
      MAGENTA.cct.p <- list(p = NA, rho1.min = NA, rho2.min = NA,rho3.min = NA)
      return(MAGENTA.cct.p)
    }
    V.diag.all <- as.matrix(Matrix::bdiag(V.diag.bydrug))
    Sigma.inv <- V.diag.all %*% U.inv %*% V.diag.all #inverse of covariance matrix of beta
  } else {
    mlist <- list()
    for (DT.ind in 1:nrow(drug_trait)) {
      d <- drug_trait[DT.ind, 2]
      k <- drug_trait[DT.ind, 1]
      mlist[[DT.ind]] <- obs.stat1[[d]][[k]]$V
    }
    U.cov.all <- V.diag.all <- Sigma.inv <- as.matrix(Matrix::bdiag(mlist))

    U.inv <- try(MASS::ginv(U.cov.all), silent = T)
    if (inherits(U.inv, "try-error")) {
      warning("The covariance matrix of U is exactly singular, MASS:ginv() function doesn't work here.")
      MAGENTA.cct.p <- list(p = NA, rho1.min = NA, rho2.min = NA, rho3.min = NA)
      return(MAGENTA.cct.p)
    }
  }

  U.all <- NULL
  for (DT.ind in 1:nrow(drug_trait)) {
    d <- drug_trait[DT.ind, 2]
    k <- drug_trait[DT.ind, 1]

    tmp <- obs.stat1[[d]][[k]]$U
    names(tmp) <- paste0(names(obs.stat1[[d]][[k]]$U), ":", k, ":",d)
    U.all <- c(U.all, tmp)
  }
  U.all <- as.matrix(U.all)

  KC <- diag(m)
  MAF <- MAF[order(match(names(MAF), snp.list))]
  if(diffweight) {
    diag.tmp <- numeric(m)
    if(length(common.SNP) != 0) {
      diag.tmp[common.SNP] <- Beta.Weights(MAF[common.SNP], weights.beta = weight.commonSNP)
      diag.tmp[-common.SNP] <- Beta.Weights(MAF[-common.SNP], weights.beta = weight.SNP)
    }else{
      diag.tmp <- Beta.Weights(MAF, weights.beta = weight.SNP)
    }

  }else{
    diag.tmp <- Beta.Weights(MAF, weights.beta = weight.SNP) #weights matrix for SNPs
  }
  diag(KC) <- diag.tmp
  WA <- diag(K)

  JC <- matrix(1, nrow = m, ncol = m)

  p.obs <- NULL
  for (ii in 1:nrow(lambda)) {
    # print(paste0("Start grid", ii, " at ", Sys.time()))
    lambdaA <- lambda[ii, 1]
    A <- WA %*% apply(KA, 2, combMTAR, rho = lambdaA) %*% WA

    lambdaC <- lambda[ii, 2]
    C <- diag(KC) %*% t(diag(KC)) * ((1 - lambdaC) * diag(m) + lambdaC * JC)
    B_D <- ar1_cor(D, lambda[ii, 3])

    B.all <- kronecker(B_D, kronecker(A, C))

    R <- chol(B.all, pivot = TRUE)
    r <- attr(R, "rank")
    if (r < nrow(B.all))
      R[(r + 1):nrow(B.all), (r + 1):nrow(B.all)] <- 0
    oo <- order(attr(R, "pivot"))
    unpivQ <- R[, oo]
    rho <- try(Get_Lambda(unpivQ %*% Sigma.inv %*% t(unpivQ)), silent = TRUE)
    if (!inherits(rho, "try-error")) {
      part1 <-t(U.all) %*% U.inv %*% V.diag.all
      Q <- as.numeric(part1 %*% B.all %*% t(part1))
      pvalues <- Get_PValue.Lambda(rho, Q)$p.value

      p.obs <- c(p.obs, pvalues)
    }
  }

  if(cct) {
    MAGENTA.cct.p <- list(p = ACAT(p.obs[which(p.obs < pval.thres)]), rho1.min = lambda[which.min(p.obs), 2],
                          rho2.min = lambda[which.min(p.obs), 1],
                          rho3.min = lambda[which.min(p.obs), 3])
  }else{
    MAGENTA.cct.p <- list(p = p.obs, rho1.min = lambda[which.min(p.obs), 2],
                          rho2.min = lambda[which.min(p.obs), 1],
                          rho3.min = lambda[which.min(p.obs), 3])
  }

  return(MAGENTA.cct.p)
}

MAGENTA.main <- function (U, V, MAF, R.C = NULL, KA = NULL, snp.list,
                          # diagTrans = FALSE,
                          zeta = NULL,
                          rho.SNP = c(0, 0.5, 1),
                          rho.trait = c(0.5, 5, 1),
                          cct = TRUE, weight.SNP = c(1, 25),
                          diffweight = FALSE,
                          threshold = 0.05,pval.thres = 1,
                          weight.commonSNP = NULL)
{
  lambda <- expand.grid(lambdaA = rho.trait, lambdaC = rho.SNP)
  K <- length(U)
  m <- length(MAF)

  if (is.null(KA)) {
    message("Without the input of genetic correlation information, an exchangeable correlation structure among traits is assumed in MAGENTA.")
    KA <- diag(K)
  }
  if (is.null(zeta)) {
    message("Without the input of zeta, MAGENTA assumes there are no overlap samples in the dataset.")
    zeta <- diag(K)
  }

  message("Conducting MAGENTA analysis ...")

  if(diffweight) {
    message("Use different weights for rare and common variants:")
    if(is.null(weight.commonSNP)){
      message(paste0("No input of weight for common variants, use the default weight of dbeta(MAF,0.5,0.5) for common variants and dbeta(MAF,", paste0(weight.SNP, collapse = ","),  ") for rare variants"))
      weight.commonSNP <- c(0.5, 0.5)
    }else{
      message(paste0("Use the weight of dbeta(MAF,", paste0(weight.SNP, collapse = ","),  ") for rare variants and the weight of dbeta(MAF,", paste0(weight.commonSNP, collapse = ","), ") for common variants"))
    }
    common.SNP <- which(MAF >= threshold)
  }

  obs.stat <- list()
  U.complete <- list()
  for (k in 1:K) {
    U.bytrait <- numeric(length(snp.list))
    V.bytrait <- matrix(0, length(snp.list), length(snp.list))
    names(U.bytrait) <- snp.list
    colnames(V.bytrait) <- rownames(V.bytrait) <- snp.list
    order1 <- order(match(names(U[[k]])[names(U[[k]]) %in%
                                          snp.list], snp.list))
    order2 <- which(snp.list %in% names(U[[k]]))
    U.bytrait[order2] <- U[[k]][order1]
    V.bytrait[order2, order2] <- V[[k]][order1, order1]
    obs.stat[[k]] <- list(U = U.bytrait, V = V.bytrait)
    U.complete[[k]] <- U.bytrait
  }

  if(is.null(R.C)) {
    R.C <- recoverLD(sumstats = list(obs.stat), snp.list = snp.list)
    message("The LD matrix is not provided and recovered from the summary statistics V.")
  }

  ## keep all non-polymorphic SNPs in the calculation, just for simplification
  obs.stat1 <- obs.stat

  if (!is.null(zeta)) {
    V.sqrt <- list()
    for (k in 1:K) {
      V.sqrt[[k]] <- diag(length(obs.stat1[[k]]$U))
      diag(V.sqrt[[k]]) <- sqrt(diag(obs.stat1[[k]]$V))
    }
    ind <- expand.grid(row.ind = 1:K, column.ind = 1:K)
    Ucov <- list()
    for (iter in 1:nrow(ind)) {
      k1 <- ind[iter, 1]
      k2 <- ind[iter, 2]

      Ucov[[iter]] <- diag(V.sqrt[[k1]])  %*%  t(diag(V.sqrt[[k2]])) * R.C

      # if(!diagTrans) {
      #   if(k1 == k2) {
      #     Ucov[[iter]] <- obs.stat1[[k1]]$V
      #   }
      # }
    }
    U.cov <- NULL
    mlist <- list()
    for (col_id in 1:K) {
      diag.ind <- which(ind$column.ind == col_id & ind$row.ind ==
                          col_id)
      U.cov.col <- NULL
      index <- which(ind$column.ind == col_id)
      row_id <- ind[index, 1]
      for (iter in 1:length(index)) {
        U.cov.col <- rbind(U.cov.col, as.matrix(zeta[row_id[iter], col_id] * Ucov[[index[iter]]]))
      }
      U.cov <- cbind(U.cov, U.cov.col)
      mlist[[col_id]] <- Ucov[[diag.ind]]
    }
    U.cov <- as.matrix(U.cov)
    U.inv <- try(MASS::ginv(U.cov), silent = T)
    if (inherits(U.inv, "try-error")) {
      warning("The covariance matrix of U is exactly singular, MASS:ginv() function doesn't work here.")
      MAGENTA.cct.p <- list(p = NA, rho1.min = NA, rho2.min = NA, det = NA, eigen = NA)
      return(MAGENTA.cct.p)
    }
    V.diag <- as.matrix(Matrix::bdiag(mlist))
    Sigma.inv <- V.diag %*% MASS::ginv(U.cov) %*% V.diag
  }else {
    mlist <- list()
    for (k in 1:K) {
      mlist[[k]] <- obs.stat1[[k]]$V
    }
    U.cov <- V.diag <- Sigma.inv <- as.matrix(Matrix::bdiag(mlist))
    U.inv <- try(MASS::ginv(U.cov), silent = T)
    if (inherits(U.inv, "try-error")) {
      warning("The covariance matrix of U is exactly singular, MASS:ginv() function doesn't work here.")
      MAGENTA.cct.p <- list(p = NA, rho1.min = NA, rho2.min = NA, det = NA, eigen = NA)
      return(MAGENTA.cct.p)
    }
  }

  U.alltraits <- NULL
  for (k in 1:K) {
    U.alltraits <- c(U.alltraits, obs.stat1[[k]]$U)
  }
  U.alltraits <- as.matrix(U.alltraits)
  KC <- diag(m)
  MAF <- MAF[order(match(names(MAF), snp.list))]
  if(diffweight) {
    diag.tmp <- numeric(m)
    if(length(common.SNP) != 0) {
      diag.tmp[common.SNP] <- Beta.Weights(MAF[common.SNP], weights.beta = weight.commonSNP)
      diag.tmp[-common.SNP] <- Beta.Weights(MAF[-common.SNP], weights.beta = weight.SNP)
    }else{
      diag.tmp <- Beta.Weights(MAF, weights.beta = weight.SNP)
    }
  }else{
    diag.tmp <- Beta.Weights(MAF, weights.beta = weight.SNP) #weights matrix for SNPs
  }

  WA <- diag(K)

  diag(KC) <- diag.tmp
  JC <- matrix(1, nrow = m, ncol = m)
  p.obs <- NULL
  for (ii in 1:nrow(lambda)) {
    lambdaA <- lambda[ii, 1]
    A <- WA %*% apply(KA, 2, combMTAR, rho = lambdaA) %*% WA

    lambdaC <- lambda[ii, 2]
    C <- diag(KC) %*% t(diag(KC)) * ((1 - lambdaC) * diag(m) + lambdaC * JC)

    B <- kronecker(A, C)
    R <- chol(B, pivot = TRUE)
    r <- attr(R, "rank")
    if (r < nrow(B))
      R[(r + 1):nrow(B), (r + 1):nrow(B)] <- 0
    oo <- order(attr(R, "pivot"))
    unpivQ <- R[, oo]
    rho <- try(Get_Lambda(unpivQ %*% Sigma.inv %*% t(unpivQ)),
               silent = TRUE)
    if (!inherits(rho, "try-error")) {
      part1 <-t(U.alltraits) %*% U.inv %*% V.diag
      Q <- as.numeric(part1 %*% B %*% t(part1))

      pvalues <- Get_PValue.Lambda(rho, Q)$p.value
      p.obs <- c(p.obs, pvalues)
    }
  }
  # p.obs <- ifelse(p.obs == 1, 0.999, p.obs)
  if(cct) {
    MAGENTA.cct.p <- list(p = ACAT(p.obs[which(p.obs < pval.thres)]),
                          rho1.min = lambda[which.min(p.obs), ]$lambdaC,
                          rho2.min = lambda[which.min(p.obs), ]$lambdaA)
  }else{
    MAGENTA.cct.p <- list(p = p.obs,
                          rho1.min = lambda[which.min(p.obs), ]$lambdaC,
                          rho2.min = lambda[which.min(p.obs), ]$lambdaA)
  }

  return(MAGENTA.cct.p)
}

MAGENTA.GEI.way1 <- function(sumstats, MAF, KA = NULL, drugStruct = NULL,
                             R.C = NULL, cct = TRUE, ref = NULL,
                             # diagTrans = FALSE,
                             diffweight = FALSE,
                             threshold = 0.05,
                             weight.commonSNP = NULL,
                             zeta = NULL, pval.thres = 1,
                             rho.SNP = c(0, 0.5, 1),
                             rho.trait = c(0.5, 5, 1),
                             rho.drug = c(0, 0.5, 1),
                             weight.SNP = c(1, 25))
{
  ## Gene-environment interaction test starting from beta, zeta is a DK by DK matrix ##
  D <- length(sumstats)
  K <- length(sumstats[[1]])
  m <- length(MAF)
  snp.list <- names(MAF)
  drug_trait <- expand.grid(trait = 1:K, drug = 1:D)

  if (is.null(drugStruct)){
    message("Assume the effects among environmental groups are independent")
    lambda <- expand.grid(lambdaA = rho.trait, lambdaC = rho.SNP, lambdaD = 0)
  }else if(drugStruct == "AR1") {
    message("Assume the first-order autoregressive structure for the effects among environmental groups")
    lambda <- expand.grid(lambdaA = rho.trait, lambdaC = rho.SNP, lambdaD = rho.drug)
  }

  if (is.null(KA)) {
    message("Without the input of genetic correlation information, an exchangeable correlation structure among traits is assumed in MAGENTA.")
    KA <- diag(K)
  }
  if (is.null(zeta)) {
    message("Without the input of zeta, MAGENTA assumes there are no overlap samples in the dataset.")
    zeta <- diag(K * D)
  }
  message("Conducting MAGENTA GEI analysis ...")

  if(diffweight) {
    message("Use different weights for rare and common variants:")
    if(is.null(weight.commonSNP)){
      message(paste0("No input of weight for common variants, use the default weight of dbeta(MAF,0.5,0.5) for common variants and dbeta(MAF,", paste0(weight.SNP, collapse = ","),  ") for rare variants"))
      weight.commonSNP <- c(0.5, 0.5)
    }else{
      message(paste0("Use the weight of dbeta(MAF,", paste0(weight.SNP, collapse = ","),  ") for rare variants and the weight of dbeta(MAF,", paste0(weight.commonSNP, collapse = ","), ") for common variants"))
    }
    common.SNP <- which(MAF >= threshold)
  }
  ## recover R.C from input of V, so MAGENTA doesn't need the extra input of R.C
  if(is.null(R.C)) {
    R.C <- recoverLD(sumstats = sumstats, snp.list = snp.list)
    message("The LD matrix is not provided and recovered from the summary statistics V.")
  }
  if(is.null(ref)) {
    ref <- D
  }
  message(paste0("Treat D = ", ref, " drug group as the reference group."))

  obs.stat <- list()
  U.complete <- list()
  for (d in 1:D) {
    U.temp <- list()
    obs.stat.temp <- list()
    for (k in 1:K) {
      U.bytrait <- numeric(length(snp.list))
      V.bytrait <- matrix(0, length(snp.list), length(snp.list))
      names(U.bytrait) <- snp.list
      colnames(V.bytrait) <- rownames(V.bytrait) <- snp.list
      order1 <- order(match(names(sumstats[[d]][[k]]$U)[names(sumstats[[d]][[k]]$U) %in% snp.list], snp.list))
      order2 <- which(snp.list %in% names(sumstats[[d]][[k]]$U))
      U.bytrait[order2] <- sumstats[[d]][[k]]$U[order1]
      V.bytrait[order2, order2] <- sumstats[[d]][[k]]$V[order1, order1]
      U.bytrait[U.bytrait == 0] <- NA
      V.bytrait[V.bytrait == 0] <- NA
      obs.stat.temp[[k]] <- list(U = U.bytrait, V = V.bytrait)
      U.temp[[k]] <- U.bytrait
    }
    obs.stat[[d]] <- obs.stat.temp
    U.complete[[d]] <- U.temp
  }

  # return(obs.stat)
  nonpolymorphic <- FALSE
  # pre-process non-polymorphic SNPs that are polymorphic at least in two groups.
  if (any(is.na(unlist(U.complete)))) {
    nonpoly.list <- list()
    analyze.snp <- list()
    for(k in 1:K) {
      ind <- NULL
      for(d in 1:D){
        ind <- rbind(ind, ifelse(is.na(obs.stat[[d]][[k]]$U), 1, 0))
      }
      nonpoly.list[[k]] <- which(apply(ind, 2, sum) > D-2)
      analyze.snp[[k]] <- colnames(ind)[-which(apply(ind, 2, sum)> D-2)]
    }

    # message("Remove non-polymorphic SNPs in at least two drug assignment groups")
    update.obs.stat <- list()
    for (d in 1:D) {
      obs.stat.temp <- list()
      for (k in 1:K) {
        ind <- nonpoly.list[[k]]
        if (length(ind) != 0) {
          U.tmp <- obs.stat[[d]][[k]]$U[-ind]
          V.tmp <- obs.stat[[d]][[k]]$V[-ind, -ind, drop = FALSE]
          obs.stat.temp[[k]] <- list(U = U.tmp, V = V.tmp)
        }
        else {
          obs.stat.temp[[k]] <- list(U = obs.stat[[d]][[k]]$U, V = obs.stat[[d]][[k]]$V)
        }
      }
      update.obs.stat[[d]] <- obs.stat.temp
    }
    update.R.C <- nonpolymorphic.fn(R.C, obs.stat, GEI = TRUE, nonpoly.list = nonpoly.list)

    nonpolymorphic <- TRUE
  }
  if (nonpolymorphic) {
    obs.stat1 <- update.obs.stat
  }else {
    obs.stat1 <- obs.stat
  }

  V.sqrt <- list() #contain the square root of all D by K covariance matrix V.
  for(d in 1:D) {
    V.sqrt.temp <- list()
    for(k in 1:K) {
      V.sqrt.temp[[k]] <- diag(length(obs.stat1[[d]][[k]]$U))
      diag(V.sqrt.temp[[k]]) <- sqrt(diag(obs.stat1[[d]][[k]]$V))
      rownames(V.sqrt.temp[[k]]) <- rownames(obs.stat1[[d]][[k]]$V)
      colnames(V.sqrt.temp[[k]]) <- colnames(obs.stat1[[d]][[k]]$V)
    }
    V.sqrt[[d]] <- V.sqrt.temp
  }

  trait.ind <- expand.grid(row.ind = 1:K, column.ind = 1:K)
  trait.drug.ind <- expand.grid(column.trait = 1:K, column.drug = 1:D, row.trait = 1:K, row.drug = 1:D)
  # Sigma.U.col <- Sigma.U <- NULL
  Ucov <- list()
  for(iter in 1:nrow(trait.drug.ind)) {
    d1 <- trait.drug.ind[iter, 4]
    d2 <- trait.drug.ind[iter, 2]
    k1 <- trait.drug.ind[iter, 3]
    k2 <- trait.drug.ind[iter, 1]
    zeta.row <- k1 + (d1 - 1) * K
    zeta.col <- k2 + (d2 - 1) * K

    if (nonpolymorphic) {
      Ucov[[iter]] <- zeta[zeta.row, zeta.col] * diag(V.sqrt[[d1]][[k1]])  %*%  t(diag(V.sqrt[[d2]][[k2]])) * update.R.C[[d1]][[which(trait.ind$row.ind == k1 & trait.ind$column.ind == k2)]]

    } else {
      Ucov[[iter]] <- zeta[zeta.row, zeta.col] * diag(V.sqrt[[d1]][[k1]])  %*%  t(diag(V.sqrt[[d2]][[k2]])) * R.C

    }

    # if(!diagTrans) {
    #   if(k1 == k2 & d1 == d2) {
    #     Ucov[[iter]] <- obs.stat1[[d1]][[k1]]$V
    #   }
    # }
    rownames(Ucov[[iter]]) <- rownames(obs.stat1[[d1]][[k1]]$V)
    colnames(Ucov[[iter]]) <- colnames(obs.stat1[[d2]][[k2]]$V)
  }

  V.inv <- list()
  for(d in 1:D) {
    V.inv.tmp <- list()
    for(k in 1:K) {
      V.subset.tmp <- obs.stat1[[d]][[k]]$V[!is.na(obs.stat1[[d]][[k]]$U),
                                            !is.na(obs.stat1[[d]][[k]]$U),drop = FALSE]
      V.inv.tmp[[k]] <- MASS::ginv(V.subset.tmp)
      rownames(V.inv.tmp[[k]]) <- rownames(V.subset.tmp)
      colnames(V.inv.tmp[[k]]) <- colnames(V.subset.tmp)

    }
    V.inv[[d]] <- V.inv.tmp
  }

  delta <- NULL
  for(d in setdiff(1:D, ref)) {
    for(k in 1:K) {
      delta.part1 <- delta.part2 <-  rep(NA, length(analyze.snp[[k]]))
      names(delta.part1) <- names(delta.part2) <- analyze.snp[[k]]

      delta.part1.tmp <- t(V.inv[[d]][[k]] %*% obs.stat1[[d]][[k]]$U[!is.na(obs.stat1[[d]][[k]]$U)])
      delta.part2.tmp <- t(V.inv[[ref]][[k]] %*% obs.stat1[[ref]][[k]]$U[!is.na(obs.stat1[[ref]][[k]]$U)])

      delta.part1[match(colnames(delta.part1.tmp), names(delta.part1))] <- delta.part1.tmp
      delta.part2[match(colnames(delta.part2.tmp), names(delta.part2))] <- delta.part2.tmp

      delta <- c(delta, delta.part1 - delta.part2)
    }
  }

  # delta <- ginv(V.diag.all) %*% U.all - kronecker(matrix(1, D), ginv(V.sum) %*% U.sum.bytrait)
  Sigma.delta <- NULL
  for(d1 in setdiff(1:D, ref)) {
    Sigma.delta.col <- NULL
    for(d2 in setdiff(1:D, ref)) {
      delta.cov <- NULL
      for(k1 in 1:K) {
        delta.cov.col <- NULL
        for(k2 in 1:K) {
          trait.drug.iter1 <- which(trait.drug.ind$row.trait == k1 & trait.drug.ind$column.trait == k2 & trait.drug.ind$row.drug == d1 & trait.drug.ind$column.drug == d2)
          trait.drug.iter4 <- which(trait.drug.ind$row.trait == k1 & trait.drug.ind$column.trait == k2 & trait.drug.ind$row.drug == ref & trait.drug.ind$column.drug == ref)

          part1 <- V.inv[[d1]][[k1]] %*%
            Ucov[[trait.drug.iter1]][!is.na(obs.stat1[[d1]][[k1]]$U),
                                     !is.na(obs.stat1[[d2]][[k2]]$U), drop = FALSE] %*%
            V.inv[[d2]][[k2]]
          part4 <- V.inv[[ref]][[k1]] %*%
            Ucov[[trait.drug.iter4]][!is.na(obs.stat1[[ref]][[k1]]$U),
                                     !is.na(obs.stat1[[ref]][[k2]]$U), drop = FALSE] %*%
            V.inv[[ref]][[k2]]
          part1 <- part1[rownames(part1) %in% analyze.snp[[k1]],
                         colnames(part1) %in% analyze.snp[[k2]],drop = FALSE]

          part4 <- part4[rownames(part4) %in% analyze.snp[[k1]],
                         colnames(part4) %in% analyze.snp[[k2]],drop = FALSE]
          part1.full <- part4.full <- matrix(NA, length(analyze.snp[[k1]]),
                                             length(analyze.snp[[k2]]))
          rownames(part1.full) <- rownames(part4.full) <- analyze.snp[[k1]]
          colnames(part1.full) <- colnames(part4.full) <- analyze.snp[[k2]]
          part1.full[match(rownames(part1), rownames(part1.full)),
                     match(colnames(part1), colnames(part1.full))] <- part1
          part4.full[match(rownames(part4), rownames(part4.full)),
                     match(colnames(part4), colnames(part4.full))] <- part4

          if(d1 == d2) {
            delta.cov.col <- cbind(delta.cov.col, part1.full + part4.full)
          }else{
            delta.cov.col <- cbind(delta.cov.col, part4.full)

          }
        }
        delta.cov <- rbind(delta.cov, delta.cov.col)
      }

      Sigma.delta.col <- cbind(Sigma.delta.col, delta.cov)
    }
    Sigma.delta <- rbind(Sigma.delta, Sigma.delta.col)
  }
  Sigma.delta <- Sigma.delta[!is.na(delta), !is.na(delta), drop = FALSE]
  Sigma.delta.inv <- try(MASS::ginv(Sigma.delta), silent = T)
  if (inherits(Sigma.delta.inv, "try-error")){
    warning("The covariance matrix of delta is exactly singular, MASS:ginv() function doesn't work here.")
    MAGENTA.cct.p <- list(p = NA, rho1.min = NA, rho2.min = NA, rho3.min = NA)
    return(MAGENTA.cct.p)
  }

  KC <- diag(m)
  MAF <- MAF[order(match(names(MAF), snp.list))]
  if(diffweight) {
    diag.tmp <- numeric(m)

    if(length(common.SNP) != 0) {
      diag.tmp[common.SNP] <- Beta.Weights(MAF[common.SNP], weights.beta = weight.commonSNP)
      diag.tmp[-common.SNP] <- Beta.Weights(MAF[-common.SNP], weights.beta = weight.SNP)
    }else{
      diag.tmp <- Beta.Weights(MAF, weights.beta = weight.SNP)
    }

  }else{
    diag.tmp <- Beta.Weights(MAF, weights.beta = weight.SNP) #weights matrix for SNPs
  }
  diag(KC) <- diag.tmp
  WA <- diag(K)

  JC <- matrix(1, nrow = m, ncol = m)

  p.obs <- NULL
  for (ii in 1:nrow(lambda)) {
    # print(ii)
    lambdaA <- lambda[ii, 1]
    A <- WA %*% apply(KA, 2, combMTAR, rho = lambdaA) %*% WA
    lambdaC <- lambda[ii, 2]
    C <- diag(KC) %*% t(diag(KC)) * ((1 - lambdaC) * diag(m) + lambdaC * JC)

    colnames(C) <- snp.list
    B_D <- ar1_cor(D-1, lambda[ii, 3])

    if (nonpolymorphic) {
      update.C <- nonpolymorphic.fn(C, obs.stat, GEI = TRUE, nonpoly.list = nonpoly.list)
      B.all <- NULL
      for(d1 in 1:(D-1)) {
        B.bydrug <- NULL
        for(d2 in 1:(D-1)) {
          B <- NULL
          for (col_id in 1:K) {
            B.col <- NULL
            index <- which(trait.ind$column.ind == col_id)
            row_id <- trait.ind[index, 1]
            for (iter in 1:length(index)) {
              B.col <- rbind(B.col, A[row_id[iter], col_id] *
                               update.C[[d]][[index[iter]]])
            }
            B <- cbind(B, B.col)
          }
          B.bydrug <- cbind(B.bydrug, B_D[d1, d2] * B)
        }
        B.all <- rbind(B.all, B.bydrug)
      }
    }else {
      B.all <- kronecker(B_D, kronecker(A, C))
    }

    B.all <- B.all[!is.na(delta), !is.na(delta), drop = FALSE]

    R <- chol(B.all, pivot = TRUE)
    r <- attr(R, "rank")
    if (r < nrow(B.all))
      R[(r + 1):nrow(B.all), (r + 1):nrow(B.all)] <- 0
    oo <- order(attr(R, "pivot"))
    unpivQ <- R[, oo]
    rho <- try(Get_Lambda(unpivQ %*% Sigma.delta.inv %*% t(unpivQ)), silent = TRUE)
    if (!inherits(rho, "try-error")) {
      part1 <- t(delta[!is.na(delta)]) %*% Sigma.delta.inv
      Q <- as.numeric(part1 %*% B.all %*% t(part1))


      pvalues <- Get_PValue.Lambda(rho, Q)$p.value
      p.obs <- c(p.obs, pvalues)
    }

  }

  # return(Q)
  if(cct) {
    MAGENTA.cct.p <- list(p = ACAT(p.obs[which(p.obs < pval.thres)]),
                          rho1.min = lambda[which.min(p.obs), ]$lambdaC,
                          rho2.min = lambda[which.min(p.obs), ]$lambdaA,
                          rho3.min = lambda[which.min(p.obs), ]$lambdaD)
  }else{
    MAGENTA.cct.p <- list(p = p.obs,
                          rho1.min = lambda[which.min(p.obs), ]$lambdaC,
                          rho2.min = lambda[which.min(p.obs), ]$lambdaA,
                          rho3.min = lambda[which.min(p.obs), ]$lambdaD)
  }

  return(MAGENTA.cct.p)
}

MAGENTA.joint.variant <- function(sumstats, MAF, KA = NULL, drugStruct = NULL,
                                  # diagTrans = FALSE,
                                  R.C = NULL, cct = TRUE,
                                  diffweight = FALSE,
                                  threshold = 0.05,
                                  weight.commonSNP = NULL,
                                  zeta = NULL, pval.thres = 1,
                                  rho.SNP = 0,
                                  rho.trait = c(0.5, 1.5, 1),
                                  rho.drug = c(0, 0.5, 1),
                                  weight.SNP = c(1, 25))
{
  # print(paste0("Start analysis at ", Sys.time()))
  D <- length(sumstats)
  K <- length(sumstats[[1]])
  m <- length(MAF)
  snp.list <- names(MAF)
  drug_trait <- expand.grid(trait = 1:K, drug = 1:D)
  if(K == 1) rho.trait <- 0
  lambda <- expand.grid(lambdaA = rho.trait, lambdaC = rho.SNP, lambdaD = rho.drug)

  obs.stat <- unlist(sumstats)
  U.stat <- obs.stat[seq(1, length(obs.stat), by = 2)]
  V.stat <- obs.stat[seq(2, length(obs.stat), by = 2)]
  V.sqrt1 <- sqrt(V.stat)

  if(K == 1) {
    U.inv <- 1/V.sqrt1^2
    U.inv[is.infinite(U.inv)] <- 0

    Sigma.inv <- diag(V.stat, ncol = D)
    WA <- 1
  }else{
    U.cov.bydrug <- list()
    for(d in 1:D) {
      U.cov.bydrug[[d]] <- zeta[[d]] * (V.sqrt1[1:K + K*(d - 1)] * matrix(1, K, K) *
                                          t(V.sqrt1[1:K + K*(d - 1)] * matrix(1, K, K)))

    }
    U.cov.all <- as.matrix(Matrix::bdiag(U.cov.bydrug))
    U.inv <- try(MASS::ginv(U.cov.all), silent = T)
    if (inherits(U.inv, "try-error")) {
      warning("The covariance matrix of U is exactly singular, MASS:ginv() function doesn't work here.")
      MAGENTA.cct.p <- list(p = NA, rho1.min = NA, rho2.min = NA,rho3.min = NA)
      return(MAGENTA.cct.p)
    }
    V.diag.all <- diag(V.stat)
    Sigma.inv <- V.diag.all %*% U.inv %*% V.diag.all #inverse of covariance matrix of beta
    WA <- diag(K)
  }

  KC <- Beta.Weights(MAF, weights.beta = weight.SNP)

  p.obs <- NULL
  for (ii in 1:nrow(lambda)) {
    C <- KC^2
    B_D <- ar1_cor(D, lambda[ii, 3])

    if(K == 1) {
      A <- 1
      B.all <- B_D * C
    }else{
      lambdaA <- lambda[ii, 1]
      A <-   WA %*% (sign(KA) * abs(KA)^lambdaA) %*% WA
      B.all <- kronecker(B_D, A * C)
    }

    R <- chol(B.all, pivot = TRUE)
    r <- attr(R, "rank")
    if (r < nrow(B.all))
      R[(r + 1):nrow(B.all), (r + 1):nrow(B.all)] <- 0
    oo <- order(attr(R, "pivot"))
    unpivQ <- R[, oo]
    tmp <- unpivQ %*% Sigma.inv %*% t(unpivQ)
    rho <- try(Get_Lambda(tmp), silent = TRUE)
    if (!inherits(rho, "try-error")) {
      if(K == 1) {
        Q <- as.numeric(U.stat %*% B.all %*% U.stat)
      }else{
        part1 <- t(U.stat) %*% U.inv %*% V.diag.all
        Q <- as.numeric(part1 %*% B.all %*% t(part1))
      }

      pvalues <- Get_PValue.Lambda(rho, Q)$p.value
      p.obs <- c(p.obs, pvalues)
    }
  }

  if(cct) {
    MAGENTA.cct.p <- list(p = ACAT(p.obs[which(p.obs < pval.thres)]),
                          rho1.min = lambda[which.min(p.obs), 2],
                          rho2.min = lambda[which.min(p.obs), 1],
                          rho3.min = lambda[which.min(p.obs), 3])
  }else{
    MAGENTA.cct.p <- list(p = p.obs, rho1.min = lambda[which.min(p.obs), 2],
                          rho2.min = lambda[which.min(p.obs), 1],
                          rho3.min = lambda[which.min(p.obs), 3])
  }

  return(MAGENTA.cct.p)
}


MTAR.main.variant <- function (U, V, MAF, R.C = NULL, KA = NULL, snp.list,
                               # diagTrans = FALSE,
                               zeta = NULL,
                               rho.SNP = 0, rho.trait = c(0.5, 1.5, 1),
                               cct = TRUE, weight.SNP = c(1, 25),
                               diffweight = FALSE,
                               threshold = 0.05,pval.thres = 1,
                               weight.commonSNP = NULL)
{
  K <- length(U)
  m <- length(MAF)
  if(K == 1) rho.trait <- 0
  lambda <- expand.grid(lambdaA = rho.trait, lambdaC = rho.SNP)

  U.stat <- unlist(U)
  V.stat <- unlist(V)
  V.sqrt1 <- sqrt(V.stat)
  if(K == 1) {
    U.inv <- 1/V.sqrt1^2
    U.inv[is.infinite(U.inv)] <- 0

    Sigma.inv <- V.stat
    WA <- 1
  }else{
    Ucov <- zeta * (V.sqrt1[1:K] * matrix(1, K, K) * t(V.sqrt1[1:K ] * matrix(1, K, K)))
    U.inv <- try(MASS::ginv(Ucov), silent = T)
    if (inherits(U.inv, "try-error")) {
      warning("The covariance matrix of U is exactly singular, MASS:ginv() function doesn't work here.")
      MAGENTA.cct.p <- list(p = NA, rho1.min = NA, rho2.min = NA,rho3.min = NA)
      return(MAGENTA.cct.p)
    }
    V.diag.all <- diag(V.stat)
    Sigma.inv <- V.diag.all %*% U.inv %*% V.diag.all #inverse of covariance matrix of beta
    WA <- diag(K)

  }

  KC <- Beta.Weights(MAF, weights.beta = weight.SNP)

  p.obs <- NULL
  for (ii in 1:nrow(lambda)) {

    # lambdaC <- lambda[ii, 2]
    C <- KC^2
    if(K == 1) {
      A <- 1
    }else{
      lambdaA <- lambda[ii, 1]
      A <-  WA %*% (sign(KA) * abs(KA)^lambdaA) %*% WA
    }
    B <- as.matrix(A * C)

    R <- chol(B, pivot = TRUE)
    r <- attr(R, "rank")
    if (r < nrow(B))
      R[(r + 1):nrow(B), (r + 1):nrow(B)] <- 0
    oo <- order(attr(R, "pivot"))
    unpivQ <- R[, oo]

    tmp <- unpivQ %*% Sigma.inv %*% t(unpivQ)
    rho <- try(Get_Lambda(tmp), silent = TRUE)

    if (!inherits(rho, "try-error")) {
      if(K == 1) {
        Q <- as.numeric(U.stat * B * U.stat)
      }else{
        part1 <-t(U.stat) %*% U.inv %*% V.diag.all
        Q <- as.numeric(part1 %*% B %*% t(part1))

      }

      pvalues <- Get_PValue.Lambda(rho, Q)$p.value
      p.obs <- c(p.obs, pvalues)
    }
  }
  # p.obs <- ifelse(p.obs == 1, 0.999, p.obs)
  if(cct) {
    MAGENTA.cct.p <- list(p = ACAT(p.obs[which(p.obs < pval.thres)]),
                          rho1.min = lambda[which.min(p.obs), ]$lambdaC,
                          rho2.min = lambda[which.min(p.obs), ]$lambdaA)
  }else{
    MAGENTA.cct.p <- list(p = p.obs,
                          rho1.min = lambda[which.min(p.obs), ]$lambdaC,
                          rho2.min = lambda[which.min(p.obs), ]$lambdaA)
  }

  return(MAGENTA.cct.p)
}


SPA_MAGENTA_diffU <- function(simdata, genotype, sumstats, zeta, zeta.list, zeta.GE,
                              selSNP, MAF, R.C, KA, MAC.thres, way, MACadj = FALSE){

  pval.main <- pval.joint <- pval.GEI <- NULL

  if(length(selSNP) != 0) {

    genotype <- genotype[, colnames(genotype) %in% selSNP, drop = FALSE]
    drug.list <- unique(simdata$covariates[, 3])
    D <- length(drug.list) # number of different drug assignments
    trait.list <- colnames(simdata$traits)
    K <- length(trait.list) # number of different drug response traits
    snp.list <- colnames(genotype)
    # m <- length(snp.list) # number of rare variants in a gene
    L <- length(grep("cov", colnames(simdata$covariates))) # number of covariates

    obs.stat.SPA <- list() # each drug has three traits summary statistics
    MAC.perSNP <- NULL # if 1 then the MAC passes the threshold, 0 otherwise
    for(d in 1:D) {
      # extract subjects with d drug assignment
      subject.ind <- which(simdata$covariates[, 3] == drug.list[d])
      subset <- list(traits = simdata$traits[subject.ind, ,drop = FALSE],
                     covariates = simdata$covariates[subject.ind, ,drop = FALSE],
                     genotype = genotype[subject.ind, ,drop = FALSE])
      if(any(names(simdata)  ==  "status")) {
        # if there is survival type data, add the failure indicators
        subset[["status"]] <- simdata$status[subject.ind, , drop = FALSE]
      }

      ## MAC count in each group
      MAC.perSNP <- rbind(MAC.perSNP, ifelse(colSums(genotype[subject.ind, ,drop = FALSE])> MAC.thres, 1, 0))
      if(way %in% c(1, 2, 3)) {
        # calculate summary statistics for different type of traits
        obs.stat.bydrug <- list()

        for(k in 1:K) {
          Y <- subset$traits[, k]
          valid.id <- which(!is.na(Y)) # extract valid value
          trait.dat <- Y[valid.id]
          status.ind <- which(colnames(subset$status) %in% paste0("Trait", k))
          if(length(unique(trait.dat)) <= 2) {
            if(length(unique(trait.dat)) < 2) {
              warning(paste0("There is only one value in ", trait.list[k], " with drug assignment value of ", drug.list[d]))
            }
            type <- "binary"
            trait.dat <- as.factor(trait.dat)
            levels(trait.dat) <- c(0, 1)[1:length(unique(trait.dat))]
            trait.dat <- as.numeric(as.character(trait.dat))
          }else if(length(status.ind) == 1){
            type <- "survival"
          }else{
            type <- "continuous"
          }

          geno.dat <- subset$genotype[valid.id, ]
          cov.dat <- model.matrix(trait.dat ~ as.matrix(subset$covariates[valid.id, 1:L, drop = FALSE]))

          if(type == "binary") {
            # use fastSPA package
            SPA.ret <- ScoreTest_SPA_wMeta(genos = t(geno.dat), pheno = trait.dat, cov = cov.dat, method = "fastSPA", output = "metaZ")

            Z <- sign(SPA.ret$p.value) * qnorm(abs(SPA.ret$p.value)/2, lower.tail = FALSE)
          }else{
            Z <- (sumstats[[d]][[k]]$U/sqrt(diag(sumstats[[d]][[k]]$V)))[selSNP]
          }
          # else if(type == "survival") {
          #   # use SPACox package
          #   surv.dat <- data.frame(Time = trait.dat, status = subset$status[valid.id, status.ind], subset$covariates[valid.id, -3])
          #   colnames(surv.dat) <- c("Time", "status", "cov1", "cov2")
          #
          #   ### using SPACox to obtain score summary statistics for each variant
          #   Phen.mtx <- cbind(ID = paste0("IID-",1:nrow(geno.dat)), surv.dat)
          #   Geno.mtx <- geno.dat
          #   rownames(Geno.mtx) <- paste0("IID-",1:nrow(geno.dat))
          #   obj.null = SPACox_Null_Model(Surv(Time, status)~cov1+cov2, data=Phen.mtx,
          #                                pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
          #   SPACox.res = SPACox(obj.null, Geno.mtx)
          #   Z <- sign(SPACox.res[, 5]) * qnorm(abs(SPACox.res[, 3])/2, lower.tail =F)
          # }

          names(Z) <- snp.list
          obs.stat.bydrug[[k]] <- Z
        }
        names(obs.stat.bydrug) <- trait.list
        obs.stat.SPA[[d]] <- obs.stat.bydrug
        # names(obs.stat.SPA) <- paste0("Drug", drug.list)
      }
      # rownames(MAC.perSNP) <- paste0("Drug", drug.list)

    }

    for(snp.id in 1:length(selSNP)) {

      U.comb <- V.comb <- list()
      for(k in 1:K) {
        if(way == 1) {
          # U equals to Z*sqrt(V)
          U.temp <- obs.stat.SPA[[1]][[k]][snp.id] *
            sqrt( sumstats[[1]][[k]]$V[selSNP[snp.id], selSNP[snp.id], drop = FALSE])
          V.temp <- sumstats[[1]][[k]]$V[selSNP[snp.id], selSNP[snp.id], drop = FALSE]
        }else if(way == 2){
          # U equals to Z
          U.temp <- obs.stat.SPA[[1]][[k]][snp.id]
          V.temp <- 1
        }else if(way == 3){
          # U equals to U, V equals to (U/Z)^2
          U.temp <- sumstats[[1]][[k]]$U[selSNP[snp.id]]
          V.temp <- (sumstats[[1]][[k]]$U[selSNP[snp.id]]/
                       obs.stat.SPA[[1]][[k]][snp.id])^2
        }else{
          # U equals to U and V equals to V, no SPA
          U.temp <- sumstats[[1]][[k]]$U[selSNP[snp.id]]
          V.temp <- sumstats[[1]][[k]]$V[selSNP[snp.id], selSNP[snp.id], drop = FALSE]
        }

        for(d in 2:D) {
          if(way == 1) {
            U.temp <- U.temp + obs.stat.SPA[[d]][[k]][snp.id] *
              sqrt(sumstats[[d]][[k]]$V[selSNP[snp.id], selSNP[snp.id], drop = FALSE])
            V.temp <- V.temp + sumstats[[d]][[k]]$V[selSNP[snp.id], selSNP[snp.id], drop = FALSE]
          }else if (way == 2){
            U.temp <- U.temp + obs.stat.SPA[[d]][[k]][snp.id]
            V.temp <- V.temp + 1
          }else if(way == 3){
            U.temp <- U.temp + sumstats[[d]][[k]]$U[selSNP[snp.id]]
            V.temp <- V.temp + (sumstats[[d]][[k]]$U[selSNP[snp.id]]/obs.stat.SPA[[d]][[k]][snp.id])^2
          }else{
            ## dont apply SPA adjust
            U.temp <- U.temp + sumstats[[d]][[k]]$U[selSNP[snp.id]]
            V.temp <- V.temp + sumstats[[d]][[k]]$V[selSNP[snp.id], selSNP[snp.id], drop = FALSE]
          }
        }

        U.comb[[k]] <- U.temp
        V.comb[[k]] <- V.temp
      }

      sumstats.perSNP <- sumstats.perSNP.joint <- list()
      for(d in 1:D) {

        sumstats.perSNP.tmp <- list()
        for(k in 1:K) {
          if(way == 1) {
            sumstats.perSNP.tmp[[k]] <- list(U = obs.stat.SPA[[d]][[k]][snp.id] *
                                               sqrt(sumstats[[d]][[k]]$V[selSNP[snp.id], selSNP[snp.id], drop = FALSE]),
                                             V = sumstats[[d]][[k]]$V[selSNP[snp.id], selSNP[snp.id], drop = FALSE])
          }else if(way == 2){
            sumstats.perSNP.tmp[[k]] <- list(U = obs.stat.SPA[[d]][[k]][snp.id],V = 1)
          }else if(way == 3){
            sumstats.perSNP.tmp[[k]] <- list(U = sumstats[[d]][[k]]$U[selSNP[snp.id]],
                                             V = (sumstats[[d]][[k]]$U[selSNP[snp.id]]/
                                                    obs.stat.SPA[[d]][[k]][snp.id])^2)
          }else{
            sumstats.perSNP.tmp[[k]] <- list(U = sumstats[[d]][[k]]$U[selSNP[snp.id]],
                                             V = sumstats[[d]][[k]]$V[selSNP[snp.id], selSNP[snp.id], drop = FALSE])
          }

        }
        sumstats.perSNP[[d]] <- sumstats.perSNP.tmp
        sumstats.perSNP.joint[[d]] <- sumstats.perSNP.tmp

        if(MACadj) {
          if(MAC.perSNP[d, snp.id] != 1) {

            # if MAC < 10, dont use the corresponding sumstats in the joint analysis
            sumstats.perSNP.joint[[d]] <- list()
            for(k in 1:K) {
              sumstats.perSNP.joint[[d]][[k]] <- list(U = 0, V = 0)
            }

          }
        }
      }

      combMTAR.p <- try(MTAR.main.variant(U = U.comb, V = V.comb, MAF = MAF[selSNP[snp.id]], R.C = 1,
                                          snp.list = names(selSNP)[snp.id],
                                          # diagTrans = TRUE,
                                          cct = TRUE, KA = KA, zeta = zeta,
                                          diffweight = FALSE, weight.SNP = NULL)$p, silent = T)
      singletrait.p <- list()
      for(k in 1:K){

        singletrait.tmp <- try(MTAR.main.variant(U = list(U.comb[[k]]), V = list(V.comb[[k]]),
                                                 MAF = MAF[selSNP[snp.id]],
                                                 R.C = 1, snp.list = names(selSNP)[snp.id],
                                                 weight.SNP = NULL,
                                                 # diagTrans = TRUE,
                                                 rho.trait = 1,
                                                 KA = matrix(1), zeta = matrix(1),
                                                 cct = TRUE, diffweight = FALSE)$p, silent = T)
        singletrait.p[[k]] <- singletrait.tmp
      }

      if(inherits(combMTAR.p, "try-error")|
         any(sapply(singletrait.p, function(x) inherits(x, "try-error")))){
        pval.main <- rbind(pval.main, rep(NA, 4))
      }else{
        pval.main <- rbind(pval.main, c(combMTAR.p, unlist(singletrait.p)))
      }

      if(sum(MAC.perSNP[, snp.id])==0) {
        pval.joint <- rbind(pval.joint, rep(NA, 4))

      }else{
        # joint test #
        combMTAR.p <- try(MTAR.PGx.joint.variant(sumstats = sumstats.perSNP.joint,
                                                 MAF = MAF[selSNP[snp.id]],
                                                 R.C = 1, weight.SNP = NULL,
                                                 KA = KA,
                                                 # diagTrans = TRUE,
                                                 zeta = zeta.list, drugStruct = "AR1",
                                                 rho.drug = c(0, 0.5, 1), cct = TRUE,
                                                 diffweight = FALSE)$p,silent = T)

        sumstats.bytrait.joint <- list()
        zeta.list1 <- list()
        for(k in 1:K) {
          sumstats.bytrait.tmp <- list()
          for(d in 1:D) {
            zeta.list1[[d]] <- matrix(1)
            sumstats.bytrait.tmp[[d]] <- list(list(U = sumstats.perSNP.joint[[d]][[k]]$U,
                                                   V = sumstats.perSNP.joint[[d]][[k]]$V))
          }
          sumstats.bytrait.joint[[k]] <- sumstats.bytrait.tmp
        }

        singletrait.p <- list()
        for(k in 1:K){
          singletrait.tmp <-
            try(MTAR.PGx.joint.variant(sumstats = sumstats.bytrait.joint[[k]],
                                       MAF = MAF[selSNP[snp.id]], KA = matrix(1),
                                       # diagTrans = TRUE,
                                       weight.SNP = NULL,
                                       R.C = 1, zeta = zeta.list1,
                                       drugStruct = "AR1",
                                       rho.drug = c(0, 0.5, 1), cct = TRUE,
                                       rho.trait = 1,
                                       diffweight = FALSE)$p, silent = T)

          singletrait.p[[k]] <- singletrait.tmp
        }
        if(inherits(combMTAR.p, "try-error")|
           any(sapply(singletrait.p, function(x) inherits(x, "try-error")))){
          pval.joint <- rbind(pval.joint, rep(NA, 4))

        }else{
          pval.joint <- rbind(pval.joint, c(combMTAR.p, unlist(singletrait.p)))

        }
      }


      if(sum(MAC.perSNP[, snp.id]) != D) {
        pval.GEI <- rbind(pval.GEI, rep(NA, 4))
      }else{
        # if SNP has MAC > MAC.thres across all the groups, conduct MAGENTA GEI test
        ## interaction effect test start ##
        combMTAR.p <-try(MTAR.PGx.GEI.diffbeta.variant(sumstats = sumstats.perSNP.joint,
                                                       MAF = MAF[selSNP[snp.id]], KA = KA,
                                                       # diagTrans = TRUE,
                                                       R.C = 1, zeta = zeta.GE, drugStruct = "AR1",
                                                       weight.SNP = NULL,
                                                       rho.drug = c(0, 0.5, 1), cct = TRUE,
                                                       diffweight = FALSE)$p,silent = T)
        # singletrait.p <- list()
        singletrait.p <- list()
        for(k in 1:K){
          singletrait.tmp <-
            try(MTAR.PGx.GEI.diffbeta.variant(sumstats = sumstats.bytrait.joint[[k]],
                                              MAF = MAF[selSNP[snp.id]], KA = matrix(1),
                                              # diagTrans = TRUE,
                                              R.C = 1,weight.SNP = NULL,
                                              zeta = zeta.GE1, drugStruct = "AR1",
                                              rho.drug = c(0, 0.5, 1), cct = TRUE,
                                              rho.trait = 1,
                                              diffweight = FALSE)$p, silent = T)

          singletrait.p[[k]] <- singletrait.tmp
        }

        if(inherits(combMTAR.p, "try-error")|
           any(sapply(singletrait.p, function(x) inherits(x, "try-error")))){
          pval.GEI <- rbind(pval.GEI, rep(NA, 4))


        }else{
          pval.GEI <- rbind(pval.GEI, c(combMTAR.p, unlist(singletrait.p)))

        }
      }

    }
  }

  main.perSNP <- cbind(pval.Multi = paste(as.vector(pval.main[, 1]), collapse = ":"),
                       pval.Single = paste(as.vector(pval.main[, 2:(K + 1)]), collapse = ":"))
  joint.perSNP <- cbind(pval.Multi = paste(as.vector(pval.joint[, 1]), collapse = ":"),
                        pval.Single = paste(as.vector(pval.joint[, 2:(K + 1)]), collapse = ":"))
  GEI.perSNP <- cbind(pval.Multi = paste(as.vector(pval.GEI[, 1]), collapse = ":"),
                      pval.Single = paste(as.vector(pval.GEI[, 2:(K + 1)]), collapse = ":"))
  pval <- list(main = main.perSNP, joint = joint.perSNP, GEI = GEI.perSNP)
  return(pval)
}


variantP3 <- function(sumstats, zeta.ret, MAF, R.C, KA, MAF.thres, test) {
  zeta <- zeta.ret$zeta.main
  zeta.list <- zeta.ret$zeta.joint
  zeta.GE <- zeta.ret$zeta.GE

  D <- length(sumstats)
  K <- length(sumstats[[1]])

  MAC10 <- which(MAF > MAF.thres)
  pval.main <- pval.joint <- pval.GEI <- NULL
  if(length(MAC10) != 0) {

    for(snp.id in 1:length(MAC10)) {
      sumstats.perSNP <- list()
      for(d in 1:D) {
        sumstats.perSNP.tmp <- list()
        for(k in 1:K) {
          sumstats.perSNP.tmp[[k]] <- list(U = sumstats[[d]][[k]]$U[MAC10[snp.id]],
                                           V = sumstats[[d]][[k]]$V[MAC10[snp.id], MAC10[snp.id], drop = FALSE])
        }
        sumstats.perSNP[[d]] <- sumstats.perSNP.tmp
      }


      if(test %in% "main") {
        # variant-level p-value for Main effect test #
        U.comb <- V.comb <- list()
        for(k in 1:K) {
          U.temp <- sumstats.perSNP[[1]][[k]]$U
          V.temp <- sumstats.perSNP[[1]][[k]]$V

          for(d in 2:D) {
            U.temp <- U.temp + sumstats.perSNP[[d]][[k]]$U
            V.temp <- V.temp + sumstats.perSNP[[d]][[k]]$V
          }

          U.comb[[k]] <- U.temp
          V.comb[[k]] <- V.temp
        }

        combMTAR.p <- try(MTAR.main.variant(U = U.comb, V = V.comb, MAF = MAF[MAC10[snp.id]], R.C = 1,
                                            snp.list = names(MAC10)[snp.id],
                                            rho.trait = c(0.5, 5, 1),
                                            # diagTrans = TRUE,
                                            cct = TRUE, KA = KA, zeta = zeta,
                                            diffweight = FALSE, weight.SNP = NULL)$p, silent = T)

        singletrait.p <- list()
        for(k in 1:K){

          singletrait.tmp <- try(MTAR.main.variant(U = list(U.comb[[k]]), V = list(V.comb[[k]]), MAF = MAF[MAC10[snp.id]],
                                                   R.C = 1, snp.list = names(MAC10)[snp.id],
                                                   rho.trait = c(0.5, 5, 1),
                                                   weight.SNP = NULL,
                                                   # diagTrans = TRUE,
                                                   KA = matrix(1), zeta = matrix(1), cct = TRUE, diffweight = FALSE)$p,
                                 silent = T)
          singletrait.p[[k]] <- singletrait.tmp
        }

        if(inherits(combMTAR.p, "try-error")|
           any(sapply(singletrait.p, function(x) inherits(x, "try-error")))){
          pval.main <- rbind(pval.main, rep(NA, 4))
        }else{
          pval.main <- rbind(pval.main, c(combMTAR.p, unlist(singletrait.p)))
        }

      }

      if(test %in% "joint") {
        # joint test #
        combMTAR.p <- try(MAGENTA.joint.variant(sumstats = sumstats.perSNP,  MAF = MAF[MAC10[snp.id]],
                                                R.C = 1, weight.SNP = NULL,
                                                KA = KA,
                                                # diagTrans = TRUE,
                                                zeta = zeta.list, drugStruct = "AR1",
                                                rho.drug = c(0, 0.5, 1), cct = TRUE,
                                                rho.trait = c(0.5, 5, 1),
                                                diffweight = FALSE)$p,silent = T)

        sumstats.bytrait <- list()
        zeta.list1 <- list()
        for(k in 1:K) {
          sumstats.bytrait.tmp <- list()
          for(d in 1:D) {
            zeta.list1[[d]] <- matrix(1)
            sumstats.bytrait.tmp[[d]] <- list(list(U = sumstats.perSNP[[d]][[k]]$U, V = sumstats.perSNP[[d]][[k]]$V))
          }
          sumstats.bytrait[[k]] <- sumstats.bytrait.tmp
        }

        singletrait.p <- list()
        for(k in 1:K){
          singletrait.tmp <-
            try(MAGENTA.joint.variant(sumstats = sumstats.bytrait[[k]], MAF = MAF[MAC10[snp.id]], KA = matrix(1),
                                      # diagTrans = TRUE,
                                      weight.SNP = NULL,
                                      rho.trait = c(0.5, 5, 1),
                                      R.C = 1, zeta = zeta.list1, drugStruct = "AR1",
                                      rho.drug = c(0, 0.5, 1), cct = TRUE,
                                      diffweight = FALSE)$p, silent = T)

          singletrait.p[[k]] <- singletrait.tmp
        }
        if(inherits(combMTAR.p, "try-error")|
           any(sapply(singletrait.p, function(x) inherits(x, "try-error")))){
          pval.joint <- rbind(pval.joint, rep(NA, 4))

        }else{
          pval.joint <- rbind(pval.joint, c(combMTAR.p, unlist(singletrait.p)))

        }
        # system.time({})
      }

      if(test %in% "GEI") {
        ## interaction effect test start ##
        combMTAR.p <-try(MAGENTA.GEI.way1(sumstats = sumstats.perSNP, MAF = MAF[MAC10[snp.id]], KA = KA,
                                          # diagTrans = TRUE,
                                          R.C = 1, zeta = zeta.GE, drugStruct = "AR1",weight.SNP = NULL,
                                          rho.drug = c(0, 0.5, 1), cct = TRUE,
                                          rho.trait = c(0.5, 5, 1),
                                          diffweight = FALSE)$p,silent = T)
        # singletrait.p <- list()
        singletrait.p <- list()
        for(k in 1:K){
          singletrait.tmp <-
            try(MAGENTA.GEI.way1(sumstats = sumstats.bytrait[[k]], MAF = MAF[MAC10[snp.id]], KA = matrix(1),
                                 # diagTrans = TRUE,
                                 R.C = 1,weight.SNP = NULL,
                                 zeta = zeta.GE1, drugStruct = "AR1",
                                 rho.drug = c(0, 0.5, 1), cct = TRUE,
                                 rho.trait = c(0.5, 5, 1),
                                 diffweight = FALSE)$p, silent = T)

          singletrait.p[[k]] <- singletrait.tmp
        }

        if(inherits(combMTAR.p, "try-error")|
           any(sapply(singletrait.p, function(x) inherits(x, "try-error")))){
          pval.GEI <- rbind(pval.GEI, rep(NA, 4))


        }else{
          pval.GEI <- rbind(pval.GEI, c(combMTAR.p, unlist(singletrait.p)))

        }
      }



    }

  }
  if(test %in% "main") {
    main.perSNP <- cbind(pval.Multi = paste(as.vector(pval.main[, 1]), collapse = ":"),
                         pval.Single = paste(as.vector(pval.main[, 2:(K + 1)]), collapse = ":"))
  }else{
    main.perSNP <- NULL
  }
  if(test %in% "joint") {
    joint.perSNP <- cbind(pval.Multi = paste(as.vector(pval.joint[, 1]), collapse = ":"),
                          pval.Single = paste(as.vector(pval.joint[, 2:(K + 1)]), collapse = ":"))
  }else{
    joint.perSNP <- NULL
  }
  if(test %in% "GEI"){
    GEI.perSNP <- cbind(pval.Multi = paste(as.vector(pval.GEI[, 1]), collapse = ":"),
                        pval.Single = paste(as.vector(pval.GEI[, 2:(K + 1)]), collapse = ":"))
  }else{
    GEI.perSNP <- NULL
  }

  pval <- list(main = main.perSNP, joint = joint.perSNP, GEI = GEI.perSNP)
  return(pval)
}

