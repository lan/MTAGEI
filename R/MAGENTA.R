MAGENTA.eq4.sumstats <- function(beta.sumstats.obj, MAF, KA, drugStruct = "AR1",
                                 trait.list, drug.list,test,ref = NULL,
                                 R.C, cct = TRUE,
                                 diffweight = FALSE,
                                 threshold = 0.05,
                                 weight.commonSNP = NULL,
                                 rho.SNP = c(0, 0.5, 1),
                                 rho.trait = c(0.5, 5, 1),
                                 rho.drug = c(0, 0.5, 1),
                                 weight.SNP = c(1, 25))
{
  if(test == "joint") {
    beta.est <- beta.sumstats.obj$beta.est
    beta.cov <- beta.sumstats.obj$beta.cov
    Sigma.beta.inv <- MASS::ginv(beta.cov)
  }else if(test == "GEI") {
    obs.stat <- beta.sumstats.obj$obs.stat
    nonpoly.list <- beta.sumstats.obj$nonpoly.list
    delta <- beta.sumstats.obj$delta
    Sigma.delta <- beta.sumstats.obj$Sigma.delta
    Sigma.delta.inv <- MASS::ginv(Sigma.delta)
  }


  ## start from beta.est and beta.cov
  D <- length(drug.list)
  K <- length(trait.list)
  m <- length(MAF)
  snp.list <- names(MAF)
  drug_trait <- expand.grid(trait = 1:K, drug = 1:D)

  if(is.null(ref)) {
    ref <- D
  }

  if (is.null(drugStruct)){
    # message("Assume the effects among environmental groups are independent")
    lambda <- expand.grid(lambdaA = rho.trait, lambdaC = rho.SNP, lambdaD = 0)
  }else if(drugStruct == "AR1") {
    # message("Assume the first-order autoregressive structure for the effects among environmental groups")
    lambda <- expand.grid(lambdaA = rho.trait, lambdaC = rho.SNP, lambdaD = rho.drug)
  }

  if (is.null(KA)) {
    message("Without the input of genetic correlation information, an exchangeable correlation structure among traits is assumed in MAGENTA.")
    KA <- diag(K)
  }
  # message("Conducting MAGENTA analysis ...")

  if(diffweight) {
    # message("Use different weights for rare and common variants:")
    if(is.null(weight.commonSNP)){
      # message(paste0("No input of weight for common variants, use the default weight of dbeta(MAF,0.5,0.5) for common variants and dbeta(MAF,", paste0(weight.SNP, collapse = ","),  ") for rare variants"))
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

    if(test == "joint") {
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
    }else if(test == "GEI"){
      colnames(C) <- snp.list
      B_D <- ar1_cor(D-1, lambda[ii, 3])
      trait.ind <- expand.grid(row.ind = 1:K, column.ind = 1:K)
      if(m == 1) {
        update.C <- C
      }else{
        update.C <- nonpolymorphic.fn(C, obs.stat, GEI = TRUE, nonpoly.list = nonpoly.list)
      }
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
              if(m > 1) {
                B.col <- rbind(B.col, A[row_id[iter], col_id] *
                                 update.C[[d1]][[index[iter]]])
              }else{
                B.col <- rbind(B.col, A[row_id[iter], col_id] * update.C)
              }

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
      rho <- try(Get_Lambda(unpivQ %*% Sigma.delta.inv %*% t(unpivQ)), silent = TRUE)
      # print(c(ii, class(rho)))
      if (!inherits(rho, "try-error")) {
        part1 <- t(delta[!is.na(delta)]) %*% Sigma.delta.inv
        Q <- as.numeric(part1 %*% B.all %*% t(part1))
        pvalues <- Get_PValue.Lambda(rho, Q)$p.value
        p.obs <- c(p.obs, pvalues)
      }
    }
  }
  if(cct) {
    p.obs <- p.obs[which(p.obs < 1)]
    if(length(p.obs) == 0){
      p.obs <- NA
    }else{
      p.obs <- ACAT(p.obs)
    }

  }
  MAGENTA.cct.p <- list(p = p.obs,
                        rhoS.min = lambda[which.min(p.obs), 2],
                        rhoT.min = lambda[which.min(p.obs), 1],
                        rhoE.min = lambda[which.min(p.obs), 3])

  return(MAGENTA.cct.p)
}

MAGENTA.eq4.4cat.sumstats <- function(beta.sumstats.obj, MAF, KA, drugStruct = "AR1",
                                      trait.list, drug.list,output = "omnibus",MAF.thres = NULL,
                                      R.C, cct = TRUE, MAC10.bygrp = NULL, test, ref = NULL,
                                      diffweight = FALSE,
                                      threshold = 0.05,
                                      weight.commonSNP = NULL,
                                      rho.SNP = c(0, 0.5, 1),
                                      rho.trait = c(0.5, 5, 1),
                                      rho.drug = c(0, 0.5, 1),
                                      weight.SNP = c(1, 25)){
  ## start from beta.est and beta.cov
  D <- length(drug.list)
  K <- length(trait.list)
  m <- length(MAF)
  snp.list <- names(MAF)
  drug_trait <- expand.grid(trait = 1:K, drug = 1:D)

  SPA <- ifelse(is.null(beta.sumstats.obj$SPA.ret), FALSE, TRUE)
  if(test == "joint") {
    ## MTMV p value
    MTMV.p <- try(MAGENTA.eq4.sumstats(beta.sumstats.obj = beta.sumstats.obj$noSPA.ret,
                                       MAF = MAF, test = test, ref = ref,
                                       KA = KA, R.C = R.C,
                                       drugStruct = drugStruct,
                                       trait.list = trait.list,
                                       drug.list = drug.list,
                                       cct = cct,
                                       diffweight = diffweight,
                                       threshold = threshold,
                                       weight.commonSNP = weight.commonSNP,
                                       rho.SNP = rho.SNP, rho.trait = rho.trait,
                                       rho.drug = rho.drug, weight.SNP = weight.SNP
    )$p,silent = T)

    ## STMV p value
    STMV.p <- rep(NA, K)
    for(k in 1:K){
      SNP.range <- grep( paste0(":", k, ":"), rownames(beta.sumstats.obj$noSPA.ret$beta.est))
      beta.sumstats.STMV <- list(beta.est = beta.sumstats.obj$noSPA.ret$beta.est[SNP.range, , drop = FALSE],
                                 beta.cov = beta.sumstats.obj$noSPA.ret$beta.cov[SNP.range, SNP.range, drop = FALSE])
      STMV.p[k] <- try(MAGENTA.eq4.sumstats(beta.sumstats.obj = beta.sumstats.STMV,
                                            MAF = MAF,
                                            KA = matrix(1),
                                            R.C = R.C,test = test, ref = ref,
                                            drugStruct = drugStruct,
                                            trait.list = trait.list[k],
                                            drug.list = drug.list,
                                            cct = cct,
                                            diffweight = diffweight,
                                            threshold = threshold,
                                            weight.commonSNP = weight.commonSNP,
                                            rho.SNP = rho.SNP, rho.trait = 1,
                                            rho.drug = rho.drug, weight.SNP = weight.SNP
      )$p, silent = T)

    }
    names(STMV.p) <- trait.list

    ## single-value p value
    if(is.null(MAF.thres) & is.null(MAC10.bygrp)) {
      message("Users should provide minor allele count per environmental group or specify a minor allele frequency threshold to decide which SNPs are included in the single-variant analysis. No single-variant analysis is conducted.")
    }else{
      # message("Conduct single-variant analysis")
      if(is.null(MAC10.bygrp)) {
        MAC10 <- MAF[MAF > MAF.thres]
        highLD <- try(caret::findCorrelation(R.C[colnames(R.C) %in% names(MAC10),
                                                 colnames(R.C) %in% names(MAC10)],
                                             cutoff = 0.98),silent = T)
        if(length(highLD)!=0 & !inherits(highLD, "try-error")) {
          selSNP <- names(MAC10)[-highLD]
        }else{
          selSNP <- names(MAC10)
        }

        MAC10.bygrp <- selSNP
      }

      MTSV.p <- rep(NA, length(MAC10.bygrp))
      STSV.p <- matrix(NA, length(MAC10.bygrp), K)
      for(SNP.id in 1:length(MAC10.bygrp)) {
        if(SPA) {
          ## use SPA-adjusted summary statistics
          ## MTSV p value
          beta.perSNP.temp <- beta.sumstats.obj$SPA.ret[[which(names(beta.sumstats.obj$SPA.ret) == MAC10.bygrp[SNP.id])]]
          beta.sumstats.MTSV <- list(beta.est = beta.perSNP.temp$beta.est,
                                     beta.cov = beta.perSNP.temp$beta.cov)
          MTSV.p[SNP.id] <- try(MAGENTA.eq4.sumstats(beta.sumstats.obj = beta.sumstats.MTSV,
                                                     MAF = MAF[names(MAF) == MAC10.bygrp[SNP.id]],
                                                     KA = KA,test = test, ref = ref,
                                                     R.C = R.C[rownames(R.C) == MAC10.bygrp[SNP.id],
                                                               colnames(R.C) == MAC10.bygrp[SNP.id], drop = FALSE],
                                                     drugStruct = drugStruct,
                                                     trait.list = trait.list,
                                                     drug.list = drug.list,
                                                     cct = cct,
                                                     diffweight = diffweight,
                                                     threshold = threshold,
                                                     weight.commonSNP = weight.commonSNP,
                                                     rho.SNP = 1, rho.trait = rho.trait,
                                                     rho.drug = rho.drug, weight.SNP = weight.SNP
          )$p, silent = T)
          ## STSV p value
          for(k in 1:K) {
            SNP.range <- grep(paste0(MAC10.bygrp[SNP.id], ":", k, ":"), rownames(beta.perSNP.temp$beta.est))
            beta.sumstats.STSV <- list(beta.est = beta.perSNP.temp$beta.est[SNP.range, , drop = FALSE],
                                       beta.cov = beta.perSNP.temp$beta.cov[SNP.range, SNP.range, drop = FALSE])
            STSV.p[SNP.id, k] <- try(MAGENTA.eq4.sumstats(beta.sumstats.obj = beta.sumstats.STSV,
                                                          MAF =  MAF[names(MAF) == MAC10.bygrp[SNP.id]],
                                                          KA = matrix(1),test = test, ref = ref,
                                                          R.C = R.C[rownames(R.C) == MAC10.bygrp[SNP.id],
                                                                    colnames(R.C) == MAC10.bygrp[SNP.id], drop = FALSE],
                                                          drugStruct = drugStruct,
                                                          trait.list = trait.list[k],
                                                          drug.list = drug.list,
                                                          cct = cct,
                                                          diffweight = diffweight,
                                                          threshold = threshold,
                                                          weight.commonSNP = weight.commonSNP,
                                                          rho.SNP = 1, rho.trait = 1,
                                                          rho.drug = rho.drug, weight.SNP = weight.SNP
            )$p, silent = T)

          }
        }else{
          ## use the original summary statistics
          ## MTSV p value
          SNP.range <- grep(paste0(MAC10.bygrp[SNP.id], ":"), rownames(beta.sumstats.obj$noSPA.ret$beta.est))
          beta.sumstats.MTSV <- list(beta.est = beta.sumstats.obj$noSPA.ret$beta.est[SNP.range, , drop = FALSE],
                                     beta.cov = beta.sumstats.obj$noSPA.ret$beta.cov[SNP.range, SNP.range, drop = FALSE])
          MTSV.p[SNP.id] <- try(MAGENTA.eq4.sumstats(beta.sumstats.obj = beta.sumstats.MTSV,
                                                     MAF = MAF[names(MAF) == MAC10.bygrp[SNP.id]],
                                                     KA = KA,test = test, ref = ref,
                                                     R.C = R.C[rownames(R.C) == MAC10.bygrp[SNP.id],
                                                               colnames(R.C) == MAC10.bygrp[SNP.id], drop = FALSE],
                                                     drugStruct = drugStruct,
                                                     trait.list = trait.list,
                                                     drug.list = drug.list,
                                                     cct = cct,
                                                     diffweight = diffweight,
                                                     threshold = threshold,
                                                     weight.commonSNP = weight.commonSNP,
                                                     rho.SNP = 1, rho.trait = rho.trait,
                                                     rho.drug = rho.drug, weight.SNP = weight.SNP
          )$p, silent = T)

          ## STSV p value
          for(k in 1:K) {
            SNP.range <- grep(paste0(MAC10.bygrp[SNP.id], ":", k, ":"), rownames(beta.sumstats.obj$noSPA.ret$beta.est))
            beta.sumstats.STSV <- list(beta.est = beta.sumstats.obj$noSPA.ret$beta.est[SNP.range, , drop = FALSE],
                                       beta.cov = beta.sumstats.obj$noSPA.ret$beta.cov[SNP.range, SNP.range, drop = FALSE])
            STSV.p[SNP.id, k] <- try(MAGENTA.eq4.sumstats(beta.sumstats.obj = beta.sumstats.STSV,
                                                          MAF =  MAF[names(MAF) == MAC10.bygrp[SNP.id]],
                                                          KA = matrix(1),test = test, ref = ref,
                                                          R.C = R.C[rownames(R.C) == MAC10.bygrp[SNP.id],
                                                                    colnames(R.C) == MAC10.bygrp[SNP.id], drop = FALSE],
                                                          drugStruct = drugStruct,
                                                          trait.list = trait.list[k],
                                                          drug.list = drug.list,
                                                          cct = cct,
                                                          diffweight = diffweight,
                                                          threshold = threshold,
                                                          weight.commonSNP = weight.commonSNP,
                                                          rho.SNP = 1, rho.trait = 1,
                                                          rho.drug = rho.drug, weight.SNP = weight.SNP
            )$p, silent = T)
          }
        }
      }
      names(MTSV.p) <- MAC10.bygrp
      rownames(STSV.p) <- MAC10.bygrp
      colnames(STSV.p) <- trait.list
    }
  }else if(test == "GEI") {
    if(is.null(ref)){
      ref <- D
    }
    trait.grid <- expand.grid(trait1 = 1:K, trait2 = 1:K)

    ## MTMV p value
    MTMV.p <- try(MAGENTA.eq4.sumstats(beta.sumstats.obj = beta.sumstats.obj$noSPA.ret,
                                       MAF = MAF, test = test, ref = ref,
                                       KA = KA, R.C = R.C,
                                       drugStruct = drugStruct,
                                       trait.list = trait.list,
                                       drug.list = drug.list,
                                       cct = cct,
                                       diffweight = diffweight,
                                       threshold = threshold,
                                       weight.commonSNP = weight.commonSNP,
                                       rho.SNP = rho.SNP, rho.trait = rho.trait,
                                       rho.drug = rho.drug, weight.SNP = weight.SNP
    )$p,silent = T)

    ## STMV p value
    STMV.p <- rep(NA, K)
    for(k in 1:K){
      SNP.range <- grep( paste0(":", k, ":"), names(beta.sumstats.obj$noSPA.ret$delta))
      beta.sumstats.STMV <- list(delta = beta.sumstats.obj$noSPA.ret$delta[SNP.range],
                                 Sigma.delta = beta.sumstats.obj$noSPA.ret$Sigma.delta[SNP.range, SNP.range, drop = FALSE],
                                 obs.stat = sapply(beta.sumstats.obj$noSPA.ret$obs.stat, function(x) x[k], simplify = FALSE),
                                 nonpoly.list = beta.sumstats.obj$noSPA.ret$nonpoly.list[k])
      STMV.p[k] <- try(MAGENTA.eq4.sumstats(beta.sumstats.obj = beta.sumstats.STMV,
                                            MAF = MAF,
                                            KA = matrix(1),
                                            R.C = R.C,test = test, ref = ref,
                                            drugStruct = drugStruct,
                                            trait.list = trait.list[k],
                                            drug.list = drug.list,
                                            cct = cct,
                                            diffweight = diffweight,
                                            threshold = threshold,
                                            weight.commonSNP = weight.commonSNP,
                                            rho.SNP = rho.SNP, rho.trait = 1,
                                            rho.drug = rho.drug, weight.SNP = weight.SNP
      )$p, silent = T)

    }
    names(STMV.p) <- trait.list
    ## single-value p value
    if(is.null(MAF.thres) & is.null(MAC10.bygrp)) {
      message("Users should provide minor allele count per environmental group or specify a minor allele frequency threshold to decide which SNPs are included in the single-variant analysis. No single-variant analysis is conducted.")
    }else{
      # message("Conduct single-variant analysis")
      if(is.null(MAC10.bygrp)) {
        MAC10 <- MAF[MAF > D * MAF.thres]
        highLD <- try(caret::findCorrelation(R.C[colnames(R.C) %in% names(MAC10),
                                                 colnames(R.C) %in% names(MAC10)],
                                             cutoff = 0.98),silent = T)
        if(length(highLD)!=0 & !inherits(highLD, "try-error")) {
          selSNP <- names(MAC10)[-highLD]
        }else{
          selSNP <- names(MAC10)
        }

        MAC10.bygrp <- selSNP
      }

      MTSV.p <- rep(NA, length(MAC10.bygrp))
      STSV.p <- matrix(NA, length(MAC10.bygrp), K)
      for(SNP.id in 1:length(MAC10.bygrp)) {
        if(SPA) {
          beta.sumstats.perSNP <- beta.sumstats.obj$SPA.ret[[which(names(beta.sumstats.obj$SPA.ret) == MAC10.bygrp[SNP.id])]]
          ## MTSV p value
          beta.sumstats.MTSV <- list(delta = beta.sumstats.perSNP$delta,
                                     Sigma.delta = beta.sumstats.perSNP$Sigma.delta,
                                     obs.stat = beta.sumstats.perSNP$obs.stat,
                                     nonpoly.list = beta.sumstats.perSNP$nonpoly.list)
          MTSV.p[SNP.id] <- try(MAGENTA.eq4.sumstats(beta.sumstats.obj = beta.sumstats.MTSV,
                                                     MAF = MAF[names(MAF) == MAC10.bygrp[SNP.id]],
                                                     KA = KA,test = test, ref = ref,
                                                     R.C = R.C[rownames(R.C) == MAC10.bygrp[SNP.id],
                                                               colnames(R.C) == MAC10.bygrp[SNP.id], drop = FALSE],
                                                     drugStruct = drugStruct,
                                                     trait.list = trait.list,
                                                     drug.list = drug.list,
                                                     cct = cct,
                                                     diffweight = diffweight,
                                                     threshold = threshold,
                                                     weight.commonSNP = weight.commonSNP,
                                                     rho.SNP = 1, rho.trait = rho.trait,
                                                     rho.drug = rho.drug, weight.SNP = weight.SNP
          )$p, silent = T)

          ## STSV p value
          for(k in 1:K) {
            SNP.range <- grep(paste0(MAC10.bygrp[SNP.id], ":", k, ":"), names(beta.sumstats.perSNP$delta))
            beta.sumstats.STSV <- list(delta = beta.sumstats.perSNP$delta[SNP.range],
                                       Sigma.delta = beta.sumstats.perSNP$Sigma.delta[SNP.range, SNP.range, drop = FALSE],
                                       obs.stat = sapply(beta.sumstats.perSNP$obs.stat, function(x) x[k], simplify = FALSE),
                                       nonpoly.list = beta.sumstats.perSNP$nonpoly.list[k])
            STSV.p[SNP.id, k] <- try(MAGENTA.eq4.sumstats(beta.sumstats.obj = beta.sumstats.STSV,
                                                          MAF =  MAF[names(MAF) == MAC10.bygrp[SNP.id]],
                                                          KA = matrix(1),test = test, ref = ref,
                                                          R.C = R.C[rownames(R.C) == MAC10.bygrp[SNP.id],
                                                                    colnames(R.C) == MAC10.bygrp[SNP.id], drop = FALSE],
                                                          drugStruct = drugStruct,
                                                          trait.list = trait.list[k],
                                                          drug.list = drug.list,
                                                          cct = cct,
                                                          diffweight = diffweight,
                                                          threshold = threshold,
                                                          weight.commonSNP = weight.commonSNP,
                                                          rho.SNP = 1, rho.trait = 1,
                                                          rho.drug = rho.drug, weight.SNP = weight.SNP
            )$p, silent = T)
          }
        }else{
          ## MTSV p value
          obs.stat.SV <- list()
          for(d in 1:D) {
            obs.stat.SV.tmp <- list()
            for(k in 1:K) {
              obs.stat.SV.tmp[[k]] <- list(U = beta.sumstats.obj$noSPA.ret$obs.stat[[d]][[k]]$U[names(beta.sumstats.obj$noSPA.ret$obs.stat[[d]][[k]]$U) == MAC10.bygrp[SNP.id]],
                                           V = beta.sumstats.obj$noSPA.ret$obs.stat[[d]][[k]]$V[names(beta.sumstats.obj$noSPA.ret$obs.stat[[d]][[k]]$U) == MAC10.bygrp[SNP.id],
                                                                                      names(beta.sumstats.obj$noSPA.ret$obs.stat[[d]][[k]]$U) == MAC10.bygrp[SNP.id], drop = FALSE])
            }
            obs.stat.SV[[d]] <- obs.stat.SV.tmp
          }
          SNP.range <- grep(paste0(MAC10.bygrp[SNP.id], ":"), names(beta.sumstats.obj$noSPA.ret$delta))
          beta.sumstats.MTSV <- list(delta = beta.sumstats.obj$noSPA.ret$delta[SNP.range],
                                     Sigma.delta = beta.sumstats.obj$noSPA.ret$Sigma.delta[SNP.range, SNP.range, drop = FALSE],
                                     obs.stat = obs.stat.SV,
                                     nonpoly.list = beta.sumstats.obj$noSPA.ret$nonpoly.list)
          MTSV.p[SNP.id] <- try(MAGENTA.eq4.sumstats(beta.sumstats.obj = beta.sumstats.MTSV,
                                                     MAF = MAF[names(MAF) == MAC10.bygrp[SNP.id]],
                                                     KA = KA,test = test, ref = ref,
                                                     R.C = R.C[rownames(R.C) == MAC10.bygrp[SNP.id],
                                                               colnames(R.C) == MAC10.bygrp[SNP.id], drop = FALSE],
                                                     drugStruct = drugStruct,
                                                     trait.list = trait.list,
                                                     drug.list = drug.list,
                                                     cct = cct,
                                                     diffweight = diffweight,
                                                     threshold = threshold,
                                                     weight.commonSNP = weight.commonSNP,
                                                     rho.SNP = 1, rho.trait = rho.trait,
                                                     rho.drug = rho.drug, weight.SNP = weight.SNP
          )$p, silent = T)

          ## STSV p value
          for(k in 1:K) {
            SNP.range <- grep(paste0(MAC10.bygrp[SNP.id], ":", k, ":"), names(beta.sumstats.obj$noSPA.ret$delta))
            beta.sumstats.STSV <- list(delta = beta.sumstats.obj$noSPA.ret$delta[SNP.range],
                                       Sigma.delta = beta.sumstats.obj$noSPA.ret$Sigma.delta[SNP.range, SNP.range, drop = FALSE],
                                       obs.stat = sapply(obs.stat.SV, function(x) x[k], simplify = FALSE),
                                       nonpoly.list = beta.sumstats.obj$noSPA.ret$nonpoly.list[k])
            STSV.p[SNP.id, k] <- try(MAGENTA.eq4.sumstats(beta.sumstats.obj = beta.sumstats.STSV,
                                                          MAF =  MAF[names(MAF) == MAC10.bygrp[SNP.id]],
                                                          KA = matrix(1),test = test, ref = ref,
                                                          R.C = R.C[rownames(R.C) == MAC10.bygrp[SNP.id],
                                                                    colnames(R.C) == MAC10.bygrp[SNP.id], drop = FALSE],
                                                          drugStruct = drugStruct,
                                                          trait.list = trait.list[k],
                                                          drug.list = drug.list,
                                                          cct = cct,
                                                          diffweight = diffweight,
                                                          threshold = threshold,
                                                          weight.commonSNP = weight.commonSNP,
                                                          rho.SNP = 1, rho.trait = 1,
                                                          rho.drug = rho.drug, weight.SNP = weight.SNP
            )$p, silent = T)
          }

        }

      }
      names(MTSV.p) <- MAC10.bygrp
      rownames(STSV.p) <- MAC10.bygrp
      colnames(STSV.p) <- trait.list
    }
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
  return(final.p)
}

MAGENTA.eq3.4cat.sumstats <- function(beta.sumstats.obj, test, MAF, KA, ref = NULL, trait.list,
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
  if(ref > D | is.null(ref)) {
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
    MTMV.p <- MAGENTA.eq3.sumstats(analyze.snp = analyze.snp[analyze.trait],
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
      STMV.p[k] <- MAGENTA.eq3.sumstats(analyze.snp = analyze.snp[k],
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
        MTSV.p[SNP.id] <- MAGENTA.eq3.sumstats(analyze.snp = lapply(analyze.snp, function(x) x[x == SNP.name])[SNP.trait],
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
            STSV.p[SNP.id, k] <- MAGENTA.eq3.sumstats(analyze.snp = lapply(analyze.snp, function(x) x[x == SNP.name])[k],
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
      message("Printing MAGENTA omnibus ", test, " p value, which combines MTMV, STMV, MTSV, STSV signal patterns")

      temp.p <- c(MTMV.p, STMV.p, MTSV.p, STSV.p)
      final.p <-  ACAT(temp.p[!is.na(temp.p) & temp.p != 1])
    }
  }else{
    message("There is no polymorphic SNP in at least 2 groups")
    final.p <- list(p = NA, MTMV = NA, STMV = NA, MTSV = NA, STSV = NA)
  }

  return(final.p)
}
MAGENTA.eq3.sumstats <- function( analyze.snp, beta.sumstats, test, MAF, KA, ref = NULL,
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
      rho <- try(Get_Lambda(unpivQ %*% Sigma.beta.inv %*% t(unpivQ)), silent = TRUE)
      if (!inherits(rho, "try-error")) {
        part1 <- t(beta) %*% Sigma.beta.inv
        Q <- as.numeric(part1 %*% B.all %*% t(part1))
        pvalues <- Get_PValue.Lambda(rho, Q)$p.value
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
      rho <- try(Get_Lambda(unpivQ %*% Sigma.delta.inv %*% t(unpivQ)), silent = TRUE)
      if (!inherits(rho, "try-error")) {
        part1 <- t(delta[!is.na(delta)]) %*% Sigma.delta.inv
        Q <- as.numeric(part1 %*% B.all %*% t(part1))
        pvalues <- Get_PValue.Lambda(rho, Q)$p.value
        p.obs <- c(p.obs, pvalues)
      }
    }
  }
  if(cct) {
    p.obs <- p.obs[which(p.obs < 1)]
    if(length(p.obs) == 0){
      p.obs <- NA
    }else{
      p.obs <- ACAT(p.obs)
    }

  }

  MTAR.cct.p <- p.obs

  return(MTAR.cct.p)
}

#' Calculate MAGENTA omnibus joint or GEI test p value using eq3/eq4
#'
#' This function allows you to calculate MAGENTA joint or GEI test p value using eq3 (start from individual-level data) or eq4 (start from summary statistics)
#' @param beta.sumstats.obj a numeric list of summary statistics generated by Get_beta_cov_UV() or Get_beta_cov_datat() functions in the MAGENTA, depending on the input data
#' @param trait.list a vector containing the traits
#' @param drug.list a vector containing the drugs or environmental groups
#' @param MAF a numeric vector containing the minor allele frequency of the \eqn{M} SNPs in the gene/SNP-set/pathway
#' @param KA a \eqn{K \times K} matrix of genetic correlation of the \eqn{K} traits
#' @param R.C a \eqn{M \times M} matrix containing LD correlation among SNPs
#' @param test if set to "joint", the function will conduct joint test; if set to "GEI", the function will return GEI p value
#' @param way should be set to "eq3" or "eq4"
#' @param MAC10.bygrp a vector containing the SNPs that users want to include in the single-variant analysis
#' @param drugStruct default is "AR1", which is auto-regressive structure with order of 1
#' @param cct whether to use Cauchy Combination test to combine the resutls, default is TRUE
#' @param MAF.thres A threshold to select SNPs with MAF > thres to be included in the single-variant analysis
#' @param ref the reference environmental group. If empty, it will be set as D (the largest group ID)
#' @param diffweight whether use different weight for common and rare variants. Default is false.
#' @param weight.SNP The weight function for each SNP. Default is dbeta(MAF, 1, 25)
#' @param weight.commonSNP if diffweight is True, the users should specify the weight function for common SNPs, like dbeta(MAF, 1, 1)
#' @param rho.SNP \eqn{rho_S} denotes the tuning parameters for between-SNP signal structure. Default is (0, 0.5,1)
#' @param rho.trait \eqn{rho_T} denotes the tuning parameters for between-trait signal structure. Default is (0.5,1, 5)
#' @param rho.drug \eqn{rho_E} denotes the tuning parameters for between-group signal structure. Default is (0, 0.5,1)
#' @param threshold a threshold to define common or rare variants. Default is 0.05.
#' @param output can be set to "omnibus", "MTMV", "STMV", "MTSV", "STSV" or "everything".
#' @return The function will return the corresponding MAGENTA omnibus p value, or multi-/single-trait multi-/single-variant p values.
#' @author Lan Luo
#' @export
#' @examples
#' \donttest{
#' #### start from individual-level data
#'#' data(rawdata)
#' names(rawdata)
#' attach(rawdata)
#'K <- 3
#'KA <- matrix(c(1, -0.03, 0.3, -0.03, 1, -0.535, 0.3, -0.535, 1), byrow = T, nrow = 3)
#'rownames(KA) <- colnames(KA) <- paste0("Trait", 1:K)
#'trait.list = paste0("Trait", 1:K)
#'drug.list = paste0("Drug", 0:1)
#'D <- 2
#'MAF <- colMeans(geno.dat[, -1])/2
#'R.C <- cor(geno.dat[, -1])
#'R.C[is.na(R.C)] <- 0
#'## start from the individual-level data ##
#'beta.sumstats.obj <- Get_beta_cov_data(geno.dat = geno.dat, dat = dat,
#' cov.list = c("cov1", "cov2"),
#' env = "treatment",
#' D = "ID", type = "continuous",
#' trait.list = paste0("Trait", 1:K),
#' SNP.list = colnames(geno.dat)[-1])
#' names(beta.sumstats.obj)
#' MAGENTA(beta.sumstats.obj = beta.sumstats.obj,
#' test = "joint", way = "eq3",MAF = MAF,
#'  R.C = R.C, KA = KA, output = "everything",
#'  trait.list = trait.list,
#'   drug.list = drug.list)
#' detach(rawdata)
#'
#' #### start from score summary statistics
#' data(sumstats.dat)
#' names(sumstats.dat)
#' attach(sumstats.dat)
#' str(common.sumstats)
#' zeta.ret <- Get_zeta(common.sumstats = common.sumstats, trait.list = paste("Trait", 1:K))
# zeta.ret
#' str(sumstats)
#' #prepare the summary statistics for joint test using trait-specific score summary statistics
#' beta.sumstats.obj <- Get_beta_cov_UV(sumstats = sumstats, MAF = MAF, R.C = R.C,
#'                                      zeta.ret = zeta.ret, KA = KA,
#'                                      test = "joint",
#'                                      trait.list = paste0("Trait", 1:K))
#' names(beta.sumstats.obj)
#' MAGENTA(beta.sumstats.obj = beta.sumstats.obj,
#'         test = "joint", way = "eq4",MAF = MAF,MAC10.bygrp = joint.MAC10,
#'         R.C = R.C, KA = KA, output = "everything",
#'         trait.list = trait.list,
#'         drug.list = drug.list)
#' detach(sumstats.dat)
#' }
MAGENTA <- function(beta.sumstats.obj, MAF, KA, drugStruct = "AR1",
                    trait.list, drug.list,output = "omnibus",MAF.thres = NULL,
                    R.C, cct = TRUE, MAC10.bygrp = NULL, test, ref = NULL,
                    diffweight = FALSE, way,
                    threshold = 0.05,
                    weight.commonSNP = NULL,
                    rho.SNP = c(0, 0.5, 1),
                    rho.trait = c(0.5, 5, 1),
                    rho.drug = c(0, 0.5, 1),
                    weight.SNP = c(1, 25)) {
  D <- length(drug.list)
  if(is.null(ref)) {
    ref <- D
  }

  if(way == "eq3") {
    message("Using eq3 to conduct MAGENTA analysis")
    final.p <- MAGENTA.eq3.4cat.sumstats(beta.sumstats.obj = beta.sumstats.obj,
                              test = test, MAF = MAF, KA = KA,
                              trait.list = trait.list,
                              output = output,
                              rho.trait = rho.trait,
                              rho.SNP = rho.SNP,
                              rho.drug = rho.drug,
                              ref = ref)
  }else if(way == "eq4") {
    message("Using eq4 to conduct MAGENTA analysis")

    final.p <- MAGENTA.eq4.4cat.sumstats(beta.sumstats.obj = beta.sumstats.obj,
                                         test = test, MAF = MAF, KA = KA, drugStruct = drugStruct,
                                         trait.list = trait.list,
                                         output = output,
                                         drug.list = drug.list,
                                         rho.trait = rho.trait,
                                         rho.SNP = rho.SNP,
                                         rho.drug = rho.drug,
                                         ref = ref,
                                         MAC10.bygrp = MAC10.bygrp,
                                         R.C = R.C)
  }
  return(final.p)
}
