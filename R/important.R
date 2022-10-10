MAGENTA.joint.sumstats <- function(beta.sumstats.obj, MAF, KA = NULL, drugStruct = "AR1",
                                   trait.list, drug.list,
                                   R.C = NULL, cct = TRUE,
                                   diffweight = FALSE,
                                   threshold = 0.05,
                                   weight.commonSNP = NULL,
                                   pval.thres = 1,
                                   rho.SNP = c(0, 0.5, 1),
                                   rho.trait = c(0.5, 5, 1),
                                   rho.drug = c(0, 0.5, 1),
                                   weight.SNP = c(1, 25))
{
  beta.est <- beta.sumstats.obj$beta.est
  beta.cov <- beta.sumstats.obj$beta.cov
  Sigma.beta.inv <- MASS::ginv(beta.cov)

  ## start from beta.est and beta.cov
  D <- length(drug.list)
  K <- length(trait.list)
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
    # rho <- try(Get_Lambda(unpivQ %*% Sigma.inv %*% t(unpivQ)), silent = TRUE)
    rho <- try(Get_Lambda(unpivQ %*% Sigma.beta.inv %*% t(unpivQ)), silent = TRUE)
    if (!inherits(rho, "try-error")) {
      # part1 <-t(U.all) %*% U.inv %*%V.diag.all
      part1 <- t(beta.est) %*% Sigma.beta.inv
      Q <- as.numeric(part1 %*% B.all %*% t(part1))
      pvalues <- Get_PValue.Lambda(rho, Q)$p.value
      # print(pvalues)
      p.obs <- c(p.obs, pvalues)
    }
  }

  if(cct) {
    MAGENTA.cct.p <- list(p = ACAT(p.obs[which(p.obs < pval.thres)]),
                          rhoS.min = lambda[which.min(p.obs), 2],
                          rhoT.min = lambda[which.min(p.obs), 1],
                          rhoE.min = lambda[which.min(p.obs), 3])
  }else{
    MAGENTA.cct.p <- list(p = p.obs,
                          rhoS.min = lambda[which.min(p.obs), 2],
                          rhoT.min = lambda[which.min(p.obs), 1],
                          rhoE.min = lambda[which.min(p.obs), 3])
  }

  return(MAGENTA.cct.p)
}

MAGENTA.joint.4cat.sumstats <- function(beta.sumstats.obj, MAF, KA = NULL, drugStruct = "AR1",
                                        trait.list, drug.list,output = "omnibus",
                                        R.C = NULL, cct = TRUE,
                                        diffweight = FALSE,
                                        threshold = 0.05,
                                        weight.commonSNP = NULL,
                                        pval.thres = 1,
                                        rho.SNP = c(0, 0.5, 1),
                                        rho.trait = c(0.5, 5, 1),
                                        rho.drug = c(0, 0.5, 1),
                                        weight.SNP = c(1, 25)){
  test <- "joint"
  # MTMV p value
  MTMV.p <- try(MAGENTA.joint.sumstats(beta.sumstats.obj = beta.sumstats.obj,
                                           MAF = MAF,
                                           KA = KA,R.C = R.C,
                                           drugStruct = drugStruct,
                                           trait.list = trait.list,
                                           drug.list = drug.list,
                                           cct = cct,
                                           diffweight = diffweight,
                                           threshold = threshold,
                                           weight.commonSNP = weight.commonSNP,
                                           pval.thres = pval.thres,
                                           rho.SNP = rho.SNP, rho.trait = rho.trait,
                                           rho.drug = rho.drug, weight.SNP = weight.SNP
  )$p,silent = T)

  STMV.p <- numeric(K)
  for(k in 1:K){
    selSNP <- grep( paste0(":", k, ":"), rownames(beta.sumstats.obj$beta.est))
    beta.sumstats.ST <- list(beta.est = beta.sumstats.obj$beta.est[selSNP, , drop = FALSE],
                             beta.cov = beta.sumstats.obj$beta.cov[selSNP, selSNP, drop = FALSE])
    STMV.p[k] <- try(MAGENTA.joint.sumstats(beta.sumstats.obj = beta.sumstats.ST,
                                                     MAF = MAF,
                                                     KA = matrix(1),
                                                     R.C = R.C,
                                                     drugStruct = drugStruct,
                                                     trait.list = trait.list[k],
                                                     drug.list = drug.list,
                                                     cct = cct,
                                                     diffweight = diffweight,
                                                     threshold = threshold,
                                                     weight.commonSNP = weight.commonSNP,
                                                     pval.thres = pval.thres,
                                                     rho.SNP = rho.SNP, rho.trait = 1,
                                                     rho.drug = rho.drug, weight.SNP = weight.SNP
    )$p, silent = T)

  }
  names(STMV.p) <- trait.list

  if(output == "everything") {
    message(paste0("Printing MAGENTA omnibus ", test, " p value as well as MTMV, STMV, MTSV, STSV ", test, " p values"))
    temp.p <- c(MTMV.p, STMV.p, MTSV.p, STSV.p)
    final.p <- list(p = ACAT(temp.p[!is.na(temp.p) & temp.p != 1]),
                    MTMV = MTMV.p,
                    STMV = STMV.p,
                    MTSV = MTSV.p,
                    STSV = STSV.p)
  }else if(output == "MTMV") {
    message("Printing MAGENTA MTMV ", test, " p value")

    final.p <- MTMV.p

  }else if(output == "STMV") {
    message("Printing MAGENTA STMV ", test, " p value")

    final.p <- STMV.p
  }else if(output == "MTSV") {
    message("Printing MAGENTA MTSV", test, " p value")

    final.p <- MTSV.p
  }else if(output == "STSV") {
    message("Printing MAGENTA STSV ", test, " p value")

    final.p <- STSV.p
  }else if(output == "omnibus") {
    message("Printing MAGENTA omnibus ", test, " p value")

    temp.p <- c(MTMV.p, STMV.p, MTSV.p, STSV.p)
    final.p <-  ACAT(temp.p[!is.na(temp.p) & temp.p != 1])
  }
return(final.p)
}

MAGENTA.calp.4cat.sumstats <- function(beta.sumstats.obj, test, MAF, KA, ref = D, trait.list,
                                       rho.trait = c(0.5, 1, 5), rho.SNP = c(0, 0.5, 1),
                                       rho.drug = c(0, 0.5, 1), output = "omnibus",
                                       weight = c(1, 25), cct = TRUE){
  beta.sumstats <- beta.sumstats.obj$beta.sumstats
  beta.sumstats.SPA <- beta.sumstats.obj$beta.sumstats.SPA
  if(test == "joint"){
    selSNP <- beta.sumstats.obj$selSNP.joint
  }else if(test == "GEI"){
    selSNP <- beta.sumstats.obj$selSNP.GEI
  }

  if(length(beta.sumstats.SPA) == 0) {
    SPA <- FALSE
  }else{
    SPA <- TRUE
  }

  D <- length(beta.sumstats)
  if(ref > D) {
    ref <- D
  }
  K <- length(trait.list)
  trait.grid <- expand.grid(trait1 = 1:K, trait2 = 1:K)
  m.total <- length(MAF)
  total.snp <- names(MAF)

  analyze.snp <- list()
  for(k in 1:K) {
    ind <- matrix(0, D, length(MAF))
    colnames(ind) <- names(MAF)
    for(d in 1:D){
      ind[d, -which(colnames(ind) %in% names(beta.sumstats[[d]]$beta.est[[k]]))] <- 1
    }
    analyze.snp[[k]] <- colnames(ind)[-which(apply(ind, 2, sum)> D-2)]
  }
  selSNP <- intersect(selSNP, Reduce(union, analyze.snp))

  beta.subset <- vector("list", D)
  for(d in 1:D) {
    beta.subset.temp <- vector("list", K)
    for(k in 1:K) {
      beta.subset.temp[[k]] <- beta.sumstats[[d]]$beta.est[[k]][names(beta.sumstats[[d]]$beta.est[[k]]) %in% analyze.snp[[k]]]
    }
    beta.subset[[d]] <- beta.subset.temp
  }
  if(length(unlist(analyze.snp)) != 0) {
    analyze.trait <- which(sapply(analyze.snp, function(x) length(x) != 0))

    ## multi-trait multi-SNP p-value
    beta.sumstats.MT <- vector("list", D)
    for(d in 1:D) {
      beta.sumstats.MT[[d]] <- list(beta.est = beta.sumstats[[d]]$beta.est[analyze.trait],
                                    beta.cov = beta.sumstats[[d]]$beta.cov[which(trait.grid$trait1 %in% analyze.trait &
                                                                                   trait.grid$trait2 %in% analyze.trait)])

    }
    MTMV.p <- MAGENTA.calp.sumstats(analyze.snp = analyze.snp[analyze.trait],
                                    beta.sumstats = beta.sumstats.MT, test = test,
                                    MAF = MAF,
                                    rho.trait = rho.trait, rho.SNP = rho.SNP,
                                    rho.drug = rho.drug,
                                    KA = KA[which(sapply(analyze.snp, function(x) length(x) != 0)),
                                            which(sapply(analyze.snp, function(x) length(x) != 0)),
                                            drop = FALSE],
                                    cct = cct, ref = ref, weight = weight)

    ## single-trait multi-SNP p-value
    STMV.p <- rep(NA, K)
    for(k in analyze.trait) {
      analyze.grp <- which(sapply(beta.subset, function(x) length(x[[k]]) != 0))
      beta.sumstats.ST <- vector("list", length(analyze.grp))
      for(d in 1:length(analyze.grp)) {
        beta.sumstats.ST[[d]] <- list(beta.est = beta.sumstats[[analyze.grp[d]]]$beta.est[k],
                                      beta.cov = beta.sumstats[[analyze.grp[d]]]$beta.cov[which(trait.grid$trait1 == k &
                                                                                                  trait.grid$trait2 == k)])
      }
      STMV.p[k] <- MAGENTA.calp.sumstats(analyze.snp = analyze.snp[k],
                                         beta.sumstats = beta.sumstats.ST,
                                         rho.trait = rho.trait, rho.SNP = rho.SNP,
                                         rho.drug = rho.drug, test = test,
                                         MAF = MAF, KA = KA[k,k, drop = FALSE],cct = cct,
                                         ref = ref, weight = weight)
    }
    names(STMV.p) <- trait.list

    if(length(selSNP) != 0) {
      MTSV.p <- numeric(length(selSNP))
      STSV.p <- matrix(NA, length(selSNP), K)
      for(SNP.id in 1:length(selSNP)) {
        SNP.name <- selSNP[SNP.id]
        SNP.MAF <- MAF[names(MAF) ==SNP.name]
        SNP.trait <- which(sapply(lapply(analyze.snp, function(x) x[x == SNP.name]),
                                  function(x) length(x) != 0))
        ## multi-trait single-SNP p-value
        beta.sumstats.SV <- vector("list", D)
        for(d in 1:D) {

          beta.sumstats.SV[[d]] <- list(beta.est = lapply(beta.sumstats[[d]]$beta.est, function(x) x[names(x) == SNP.name])[SNP.trait],
                                        beta.cov = lapply(beta.sumstats[[d]]$beta.cov, function(x) x[rownames(x) == SNP.name,
                                                                                                     colnames(x) == SNP.name,
                                                                                                     drop = FALSE])[which(trait.grid$trait1 %in% SNP.trait &
                                                                                                                            trait.grid$trait2 %in% SNP.trait)])

          if(SPA) {
            SPA.act.all <- any(sapply(beta.sumstats.SPA, function(x) sapply(x$SPA.act, function(y) y[names(y) == SNP.name])))

            if(SPA.act.all) {
              beta.cov.SPA <- SPA.convertV.bygrp(SNP.name,d, beta.sumstats = beta.sumstats,
                                                 beta.sumstats.SPA = beta.sumstats.SPA, K = K)
              beta.sumstats.SV[[d]] <- list(beta.est = lapply(beta.sumstats.SPA[[d]]$beta.est, function(x) x[names(x) == SNP.name])[SNP.trait],
                                            beta.cov = beta.cov.SPA)
            }

          }
        }
        SNP.grp <- which(sapply(beta.sumstats.SV, function(x) length(unlist(x$beta.est)) != 0)) # this SNP may be nonpolymorphic in one grp
        MTSV.p[SNP.id] <- MAGENTA.calp.sumstats(analyze.snp = lapply(analyze.snp, function(x) x[x == SNP.name])[SNP.trait],
                                                beta.sumstats = beta.sumstats.SV[SNP.grp],
                                                rho.trait = rho.trait, rho.SNP = rho.SNP,
                                                rho.drug = rho.drug, test = test,
                                                MAF = SNP.MAF, KA = KA[SNP.trait, SNP.trait, drop = FALSE],
                                                cct = cct, ref = ref, weight = weight)

        ## single-trait single-SNP p-value
        for(k in SNP.trait) {

          beta.sumstats.SVST <- vector("list", D)
          for(d in 1:D) {

            beta.sumstats.SVST[[d]] <- list(beta.est = lapply(beta.sumstats[[d]]$beta.est,
                                                              function(x) x[names(x) == SNP.name])[k],
                                            beta.cov = lapply(beta.sumstats[[d]]$beta.cov,
                                                              function(x) x[rownames(x) == SNP.name,
                                                                            colnames(x) == SNP.name,
                                                                            drop = FALSE])[which(trait.grid$trait1 == k &
                                                                                                   trait.grid$trait2 == k)])

            if(SPA){
              if(SPA.act.all) {
                beta.cov.SPA <- SPA.convertV.bygrp(SNP.name,d, beta.sumstats = beta.sumstats,
                                                   beta.sumstats.SPA = beta.sumstats.SPA, K = K)
                beta.sumstats.SVST[[d]] <- list(beta.est = lapply(beta.sumstats.SPA[[d]]$beta.est,
                                                                  function(x) x[names(x) == SNP.name])[k],
                                                beta.cov = lapply(beta.cov.SPA,
                                                                  function(x) x[rownames(x) == SNP.name,
                                                                                colnames(x) == SNP.name,
                                                                                drop = FALSE])[which(trait.grid$trait1 == k &
                                                                                                       trait.grid$trait2 == k)]
                )
              }
            }
          }
          SNP.grp.pertrait <- which(sapply(beta.sumstats.SVST, function(x) length(x$beta.est[[1]]) != 0)) # this SNP may be nonpolymorphic in one grp
          if(length(SNP.grp.pertrait) <= 1) {
            next
          }else{
            STSV.p[SNP.id, k] <- MAGENTA.calp.sumstats(analyze.snp = lapply(analyze.snp, function(x) x[x == SNP.name])[k],
                                                       beta.sumstats = beta.sumstats.SVST[SNP.grp.pertrait],
                                                       rho.trait = rho.trait, rho.SNP = rho.SNP,
                                                       rho.drug = rho.drug, test = test,
                                                       MAF = SNP.MAF, KA = KA[k,k,drop = FALSE], cct = cct,
                                                       ref = ref, weight = weight)
          }

        }
      }
      names(MTSV.p) <- rownames(STSV.p) <- selSNP
      colnames(STSV.p) <- trait.list
    }else{
      MTSV.p <- STSV.p <- NA
    }
    if(output == "everything") {
      message(paste0("Printing MAGENTA omnibus ", test, " p value as well as MTMV, STMV, MTSV, STSV ", test, " p values"))
      temp.p <- c(MTMV.p, STMV.p, MTSV.p, STSV.p)
      final.p <- list(p = ACAT(temp.p[!is.na(temp.p) & temp.p != 1]),
                      MTMV = MTMV.p,
                      STMV = STMV.p,
                      MTSV = MTSV.p,
                      STSV = STSV.p)
    }else if(output == "MTMV") {
      message("Printing MAGENTA MTMV ", test, " p value")

      final.p <- MTMV.p

    }else if(output == "STMV") {
      message("Printing MAGENTA STMV ", test, " p value")

      final.p <- STMV.p
    }else if(output == "MTSV") {
      message("Printing MAGENTA MTSV", test, " p value")

      final.p <- MTSV.p
    }else if(output == "STSV") {
      message("Printing MAGENTA STSV ", test, " p value")

      final.p <- STSV.p
    }else if(output == "omnibus") {
      message("Printing MAGENTA omnibus ", test, " p value")

      temp.p <- c(MTMV.p, STMV.p, MTSV.p, STSV.p)
      final.p <-  ACAT(temp.p[!is.na(temp.p) & temp.p != 1])
    }
  }else{
    message("There is no polymorphic SNP in at least 2 groups")
    final.p <- list(p = NA, MTMV = NA, STMV = NA, MTSV = NA, STSV = NA)
  }

  return(final.p)
}
MAGENTA.calp.sumstats <- function( analyze.snp, beta.sumstats, test, MAF, KA, ref = D,
                                   rho.trait = c(0.5, 1, 5), rho.SNP = c(0, 0.5, 1),
                                   rho.drug = c(0, 0.5, 1),
                                   weight = c(1, 25), cct = TRUE){
  # beta.sumstats <- beta.sumstats.obj$beta.sumstats
  D <- length(beta.sumstats)
  if(ref > D) {
    ref <- D
  }
  K <- length(beta.sumstats[[1]]$beta.est)
  trait.grid <- expand.grid(trait1 = 1:K, trait2 = 1:K)
  m.total <- length(MAF)
  total.snp <- names(MAF)

  if(test == "joint") {
    #### concatenate beta
    beta <- NULL
    for(d in 1:D) {
      for(k in 1:K) {
        beta.part1 <- rep(NA, length(analyze.snp[[k]]))
        names(beta.part1) <- analyze.snp[[k]]
        beta.part1.temp <- beta.sumstats[[d]]$beta.est[[k]][names(beta.sumstats[[d]]$beta.est[[k]]) %in% analyze.snp[[k]]]
        beta.part1[match(names(beta.part1.temp), names(beta.part1))] <- beta.part1.temp

        beta <- c(beta, beta.part1)
      }
    }

    trait.grid <- expand.grid(trait1 = 1:K, trait2 = 1:K)
    Sigma.beta <- list()
    for(d1 in 1:D) {
      beta.cov.bydrug <- NULL
      for(k1 in 1:K) {
        beta.cov.bydrug.row <- NULL
        for(k2 in 1:K) {
          part1 <- beta.sumstats[[d1]]$beta.cov[[which(trait.grid$trait1 == k1 & trait.grid$trait2 == k2)]]
          part1 <- part1[rownames(part1) %in% analyze.snp[[k1]],
                         colnames(part1) %in% analyze.snp[[k2]],drop = FALSE]
          part1.full <- matrix(NA, length(analyze.snp[[k1]]),
                               length(analyze.snp[[k2]]))
          rownames(part1.full) <- analyze.snp[[k1]]
          colnames(part1.full) <- analyze.snp[[k2]]
          part1.full[match(rownames(part1), rownames(part1.full)),
                     match(colnames(part1), colnames(part1.full))] <- part1
          beta.cov.bydrug.row <- cbind(beta.cov.bydrug.row, part1)
        }
        beta.cov.bydrug <- rbind(beta.cov.bydrug, beta.cov.bydrug.row)
      }
      Sigma.beta[[d1]] <- beta.cov.bydrug
    }
    Sigma.beta <- as.matrix(Matrix::bdiag(Sigma.beta))
    Sigma.beta.inv <- try(MASS::ginv(Sigma.beta), silent = T)

  }else if(test == "GEI") {
    #### calculate delta
    delta <- NULL
    for(d in setdiff(1:D, ref)) {
      for(k in 1:K) {
        beta.part1 <- beta.part2 <- rep(NA, length(analyze.snp[[k]]))
        names(beta.part1) <- names(beta.part2) <- analyze.snp[[k]]
        beta.part1.temp <- beta.sumstats[[d]]$beta.est[[k]][names(beta.sumstats[[d]]$beta.est[[k]]) %in% analyze.snp[[k]]]
        beta.part1[match(names(beta.part1.temp), names(beta.part1))] <- beta.part1.temp
        beta.part2.temp <-  beta.sumstats[[ref]]$beta.est[[k]][names(beta.sumstats[[ref]]$beta.est[[k]]) %in% analyze.snp[[k]]]
        beta.part2[match(names(beta.part2.temp), names(beta.part2))] <- beta.part2.temp

        delta <- c(delta, beta.part1 - beta.part2)
      }
    }
    # delta[is.na(delta)] <- 0

    Sigma.delta <- NULL
    for(d1 in setdiff(1:D, ref)) {
      Sigma.delta.col <- NULL
      for(d2 in setdiff(1:D, ref)) {
        delta.cov <- NULL
        for(k1 in 1:K) {
          delta.cov.col <- NULL
          for(k2 in 1:K) {
            row.col.ind <- which(trait.grid[, 1] == k1 & trait.grid[, 2] == k2)
            part1 <- beta.sumstats[[d1]]$beta.cov[[row.col.ind]]
            part4 <- beta.sumstats[[ref]]$beta.cov[[row.col.ind]]

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
            # print(delta.cov.col)
          }
          delta.cov <- rbind(delta.cov, delta.cov.col)
        }

        Sigma.delta.col <- cbind(Sigma.delta.col, delta.cov)
      }
      Sigma.delta <- rbind(Sigma.delta, Sigma.delta.col)
    }
    Sigma.delta <- Sigma.delta[!is.na(delta), !is.na(delta), drop = FALSE]
    Sigma.delta.inv <- try(MASS::ginv(Sigma.delta), silent = T)

  }

  trait.ind <- expand.grid(row.ind = 1:K, column.ind = 1:K)
  trait.drug.ind <- expand.grid(column.trait = 1:K, column.drug = 1:D, row.trait = 1:K, row.drug = 1:D)
  lambda <- expand.grid(rho.trait =rho.trait, rho.SNP = rho.SNP, rho.drug = rho.drug)

  KC <- diag(m.total)
  diag.tmp <- Beta.Weights(MAF, weights.beta = weight) #weights matrix for SNPs

  diag(KC) <- diag.tmp
  WA <- diag(K)
  JC <- matrix(1, nrow = m.total, ncol = m.total)

  if(test == "joint") {
    p.obs <- NULL
    for (ii in 1:nrow(lambda)) {
      rho.trait <- lambda[ii, 1]
      A <- WA %*% apply(KA, 2, combMTAR, rho = rho.trait) %*% WA
      rho.SNP <- lambda[ii, 2]
      C <- diag(KC) %*% t(diag(KC)) * ((1 - rho.SNP) * diag(m.total) + rho.SNP * JC)
      colnames(C) <- rownames(C) <- total.snp
      B_D <- ar1_cor(D, lambda[ii, 3])

      B.all <- NULL
      for(d1 in 1:D) {
        B.bydrug <- NULL
        for(d2 in 1:D) {
          B <- NULL
          for (col_id in 1:K) {
            B.col <- NULL
            index <- which(trait.ind$column.ind == col_id)
            row_id <- trait.ind[index, 1]
            for (iter in 1:length(index)) {
              C1 <-  remove.nonpoly(C, nonpoly.k1 = setdiff(total.snp, analyze.snp[[row_id[iter]]]),
                                    nonpoly.k2 = setdiff(total.snp, analyze.snp[[col_id]]))
              B.col <- rbind(B.col, A[row_id[iter], col_id] * C1 ## remove nonpolymorphic SNPs
              )
            }
            B <- cbind(B, B.col)
          }
          B.bydrug <- cbind(B.bydrug, B_D[d1, d2] * B)
        }
        B.all <- rbind(B.all, B.bydrug)
      }
      # B.all <- B.all[!is.na(beta), !is.na(beta), drop = FALSE]
      R <- chol(B.all, pivot = TRUE)
      r <- attr(R, "rank")
      if (r < nrow(B.all))
        R[(r + 1):nrow(B.all), (r + 1):nrow(B.all)] <- 0
      oo <- order(attr(R, "pivot"))
      unpivQ <- R[, oo]
      rho <- try(SKAT:::Get_Lambda(unpivQ %*% Sigma.beta.inv %*% t(unpivQ)), silent = TRUE)
      if (!inherits(rho, "try-error")) {
        part1 <- t(beta) %*% Sigma.beta.inv
        Q <- as.numeric(part1 %*% B.all %*% t(part1))
        pvalues <- SKAT:::Get_PValue.Lambda(rho, Q)$p.value
        p.obs <- c(p.obs, pvalues)
      }
    }
  }else if (test == "GEI") {
    p.obs <- NULL
    for (ii in 1:nrow(lambda)) {
      rho.trait <- lambda[ii, 1]
      A <- WA %*% apply(KA, 2, combMTAR, rho = rho.trait) %*% WA
      rho.SNP <- lambda[ii, 2]
      C <- diag(KC) %*% t(diag(KC)) * ((1 - rho.SNP) * diag(m.total) + rho.SNP * JC)
      colnames(C) <- rownames(C) <- total.snp
      B_D <- ar1_cor(D-1, lambda[ii, 3])

      B.all <- NULL
      for(d1 in setdiff(1:D, ref)) {
        B.bydrug <- NULL
        for(d2 in setdiff(1:D, ref)) {
          B <- NULL
          for (col_id in 1:K) {
            B.col <- NULL
            index <- which(trait.ind$column.ind == col_id)
            row_id <- trait.ind[index, 1]
            for (iter in 1:length(index)) {
              C1 <-  remove.nonpoly(C, nonpoly.k1 = setdiff(total.snp, analyze.snp[[row_id[iter]]]),
                                    nonpoly.k2 = setdiff(total.snp, analyze.snp[[col_id]]))
              B.col <- rbind(B.col, A[row_id[iter], col_id] * C1 ## remove nonpolymorphic SNPs
              )
            }
            B <- cbind(B, B.col)
          }
          B.bydrug <- cbind(B.bydrug, B_D[d1, d2] * B)
        }
        B.all <- rbind(B.all, B.bydrug)
      }
      B.all <- B.all[!is.na(delta), !is.na(delta), drop = FALSE]
      R <- chol(B.all, pivot = TRUE)
      r <- attr(R, "rank")
      if (r < nrow(B.all))
        R[(r + 1):nrow(B.all), (r + 1):nrow(B.all)] <- 0
      oo <- order(attr(R, "pivot"))
      unpivQ <- R[, oo]
      rho <- try(SKAT:::Get_Lambda(unpivQ %*% Sigma.delta.inv %*% t(unpivQ)), silent = TRUE)
      if (!inherits(rho, "try-error")) {
        part1 <- t(delta[!is.na(delta)]) %*% Sigma.delta.inv
        Q <- as.numeric(part1 %*% B.all %*% t(part1))
        pvalues <- SKAT:::Get_PValue.Lambda(rho, Q)$p.value
        p.obs <- c(p.obs, pvalues)
      }
    }
  }

  if(cct) {
    MTAR.cct.p <- ACAT(p.obs[which(p.obs < 1)])
  }else{
    MTAR.cct.p <- p.obs
  }
  return(MTAR.cct.p)
}
