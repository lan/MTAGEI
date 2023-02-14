cal.sum <- function(simdata, genotype, covariance = TRUE){

  drug.list <- unique(simdata$covariates[, 3])
  D <- length(drug.list) # number of different drug assignments
  trait.list <- colnames(simdata$traits)
  K <- length(trait.list) # number of different drug response traits
  snp.list <- colnames(genotype)
  m <- length(snp.list) # number of rare variants in a gene
  L <- length(grep("cov", colnames(simdata$covariates))) # number of covariates

  obs.stat <- list() # each drug has three traits summary statistics
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

      geno.dat <- subset$genotype[valid.id, ,drop = FALSE]
      cov.dat <- model.matrix(trait.dat ~ as.matrix(subset$covariates[valid.id, 1:L, drop = FALSE]))
      if(type == "continuous") {
        # message(paste0("Calculating summary statistics for continuous type Trait", k, " with drug assignment value of ", drug.list[d]))

        mod <- summary(lm(trait.dat ~ cov.dat))
        s2 <- mod$sigma^2
        resid <- mod$residuals
        U <- as.vector(crossprod(geno.dat, resid)) / s2 #as.vector(t(geno.dat) %*% resid) / s2
        names(U) <- snp.list
        if(covariance) {
          W <- crossprod(geno.dat, geno.dat) - crossprod(geno.dat, cov.dat) %*%
            solve(crossprod(cov.dat, cov.dat)) %*% crossprod(cov.dat, geno.dat)
          V <- W / s2
          if(m > 2) {
            colnames(V) <- rownames(V) <- snp.list
          }else{
            names(V) <- snp.list
          }
        }else{
          V <- numeric(m)
          for(j in 1:m){
            G <- geno.dat[, j]

            V[j] <- (crossprod(G, G) - crossprod(G, cov.dat) %*%
                       solve(crossprod(cov.dat, cov.dat)) %*% crossprod(cov.dat, G))/s2
            # V[j] <- as.numeric(t(G) %*% G - (t(G) %*% cov.dat) %*%
            #                      solve(t(cov.dat) %*% cov.dat) %*%
            #                      (t(cov.dat) %*% G))/s2
          }
          names(V) <- snp.list
        }
      }else if(type == "binary"){

        # message(paste0("Calculating summary statistics for binary type Trait", k, " with drug assignment value of ", drug.list[d]))
        mod <- glm(trait.dat ~ cov.dat, family = "binomial")
        b1 <- mod$fitted.values
        b2 <- b1 * (1 - b1)
        b2.prod <- b2 %*% matrix(1, 1, ncol(cov.dat)) * cov.dat

        U <- crossprod(geno.dat, (trait.dat - b1))
        U <- as.vector(U)
        names(U) <- snp.list

        if(covariance) {
          Sigma.inv <- diag(as.vector(b2), length(b2))
          Sigma.cov <- crossprod(Sigma.inv, cov.dat)

          V.tmp <- Sigma.inv - Sigma.cov %*% tcrossprod(solve(crossprod(cov.dat, Sigma.cov)), Sigma.cov)
          V <- crossprod(geno.dat, V.tmp) %*% geno.dat

          if(m > 2) {
            colnames(V) <- rownames(V) <- snp.list
          }else{
            names(V) <- snp.list
          }
        }else{

          V <- numeric(m)
          for(j in 1:m){
            G <- as.matrix(geno.dat[, j])
            V[j] <- crossprod(G, (b2 * G)) - crossprod((b2 * G), cov.dat) %*%
              solve(crossprod(cov.dat, b2.prod)) %*%
              crossprod(b2.prod, G)
          }

          names(U) <- names(V) <- snp.list
        }
      }
      obs.stat.bydrug[[k]] <- list(U = U, V = V)
    }
    names(obs.stat.bydrug) <- trait.list
    obs.stat[[d]] <- obs.stat.bydrug
  }
  names(obs.stat) <- paste0("Drug", drug.list)
  return(obs.stat)
}
#' Calculate the summary statistics required by MTAGEI given the score summary statistics
#'
#' This function allows you to calculate the summary statistics (beta, beta.cov) given the score summary statistics.
#' @param sumstats.obj a list of score summary statistics (U, V) in each environmental group for each trait, if the traits are binary, should also include the SPA-adjusted score summary statistics
#' @param MAF a numeric vector containing the minor allele frequency of the \eqn{M} SNPs in the gene/SNP-set/pathway
#' @param KA a \eqn{K \times K} matrix of genetic correlation of the \eqn{K} traits
#' @param MAF.thres A threshold to select SNPs with MAF > thres to be included in the single-variant analysis
#' @param R.C a \eqn{M \times M} matrix containing LD correlation among SNPs
#' @param trait.list a vector containing the traits
#' @param zeta.ret a list of zeta matrix to correct for overlapping samples, returned by function Get_zeta
#' @param ref the reference environmental group. If empty, it will be set as D (the largest group ID)
#' @param main.thres a threshold to prompt warning if the genetic main effect is too strong. Default is 5e-3
#' @param test if set to "joint", the function will conduct joint test; if set to "GEI", the function will return GEI p value
#' @param type should be set to "continuous" or "binary" for the traits type.
#' @return summary statistics (beta, beta.cov) required by MTAGEI
#' @author Lan Luo
#' @export
#' @examples
#' \donttest{
#' data("sumstats.dat")
#' names(sumstats.dat)
#' attach(sumstats.dat)
#' str(common.sumstats)
#' K <- 3
#' zeta.ret <- Get_zeta(common.sumstats = common.sumstats, trait.list = paste("Trait", 1:K))
#' # zeta.ret
#' str(sumstats)
#' ##prepare the summary statistics for joint test using trait-specific score summary statistics
#' beta.sumstats.obj <- Get_beta_cov_UV(sumstats.obj = sumstats.dat, MAF = MAF, R.C = R.C,
#'                                      zeta.ret = zeta.ret, KA = KA,
#'                                      test = "joint",  type = "continuous",
#'                                    trait.list = paste0("Trait", 1:K))
#' names(beta.sumstats.obj)
#' detach(sumstats.dat)
#' }
Get_beta_cov_UV <- function(sumstats.obj, MAF, zeta.ret, R.C, KA, MAF.thres = 0.000588,
                            trait.list, ref = NULL, type,
                            main.thres = 5e-3,
                            test){

  sumstats <- sumstats.obj$sumstats
  SPA <- ifelse(type == "binary", TRUE, FALSE)
  if(SPA){
    sumstats.SPA <- sumstats.obj$sumstats.SPA
  }
  D <- length(sumstats)
  K <- length(sumstats[[1]])
  m <- length(MAF)
  snp.list <- names(MAF)
  drug_trait <- expand.grid(trait = 1:K, drug = 1:D)

  if(is.null(ref)){
    ref <- D
  }

  # recover R.C from input of V, so MTAGEI doesn't need the extra input of R.C
  if(is.null(R.C)) {
    R.C <- recoverLD(sumstats = sumstats, snp.list = snp.list)
    message("The LD matrix is not provided and recovered from the summary statistics V.")
  }
  ## first calculate main p value
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
  combMTAR.p <- try(MTAGEI.main(U = U.comb, V = V.comb, MAF = MAF, R.C = R.C, snp.list = names(MAF),
                                 cct = TRUE, KA = KA, zeta = zeta.ret$zeta.main, diffweight = FALSE)$p,
                    silent = T)
  singletrait.p <- list()
  for(k in 1:K){
    singletrait.p[[k]] <- try(MTAGEI.main(U = list(U.comb[[k]]), V = list(V.comb[[k]]), MAF = MAF, R.C = R.C,
                                           snp.list = names(MAF), KA = matrix(1),
                                           rho.trait = c(0.5, 5, 1),
                                           zeta = matrix(1), cct = TRUE, diffweight = FALSE)$p,silent = T)


  }
  if(inherits(combMTAR.p, "try-error")|
     any(sapply(singletrait.p, function(x) inherits(x, "try-error")))){
    main.p <- rep(NA, 4)
  }else{
    main.p <- c(combMTAR.p, unlist(singletrait.p))
  }
  names(main.p) <- c("Multi.M", paste0(trait.list, ".M"))
  # print(main.p)
  if(SPA) {
    variant.p <- variantP3(sumstats = sumstats.SPA, zeta.ret = zeta.ret,
                           MAF = MAF, R.C = R.C, KA = KA, MAF.thres = MAF.thres,
                           test = "main", SPA = TRUE)
  }else{
    variant.p <- variantP3(sumstats = sumstats, zeta.ret = zeta.ret,
                           MAF = MAF, R.C = R.C, KA = KA, MAF.thres = MAF.thres,
                           test = "main", SPA = FALSE)
  }

  p.main.wo <- main.p
  main.varp <- cbind(as.numeric(stringr::str_split(variant.p$main[1], ":")[[1]]),
                     matrix(as.numeric(stringr::str_split(variant.p$main[2], ":")[[1]]), ncol = 3))
  MAC10 <- which(MAF > MAF.thres)
  rownames(main.varp) <- names(MAC10)
  highLD <- try(caret::findCorrelation(R.C[colnames(R.C) %in% names(MAC10), colnames(R.C) %in% names(MAC10)],
                                       cutoff = 0.98),silent = T)
  if(length(highLD)!=0 & !inherits(highLD, "try-error")) {
    main.new <- main.varp[-highLD, ,drop = FALSE]
  }else{
    main.new <- main.varp
  }

  p.main.w <- ACAT(c(p.main.wo, main.new))
  if(p.main.w < main.thres) {
    warning(paste0("The genetic main effect is strong (p-value < ", main.thres,
                   "), the GEI p-value only using score summary statistics may be biased. Users should consider provide individual-level data to obtain a more accurate estimates for GEI analysis"))
  }
  noSPA.ret <- convert_UV_to_betacov(D = D, K = K, snp.list = snp.list, sumstats = sumstats, R.C = R.C, test = test, ref = ref)
  if(SPA){
    SPA.ret <- list()
    for(snp.id in 1:length(snp.list)) {
      sumstats.perSNP <- list()
      for(d in 1:D) {
        sumstats.perSNP.temp <- list()
        for(k in 1:K) {
          sumstats.perSNP.temp[[k]] <- list(U = sumstats.SPA[[d]][[k]]$U[snp.id],
                                            V = matrix(sumstats.SPA[[d]][[k]]$V[snp.id]))
        }
        sumstats.perSNP[[d]] <- sumstats.perSNP.temp
      }
      SPA.ret[[snp.id]] <- convert_UV_to_betacov(D = D, K = K, snp.list = snp.list[snp.id],
                                                 sumstats = sumstats.perSNP, R.C = R.C[snp.id, snp.id, drop = FALSE], test = test, ref = ref)
    }
    names(SPA.ret) <- snp.list
  }else{
    SPA.ret <- NULL
  }
  final.ret <- list(noSPA.ret = noSPA.ret, SPA.ret = SPA.ret)
  return(final.ret)
}

convert_UV_to_betacov <- function(D, K, snp.list, sumstats, test, R.C, ref){
  drug_trait <- expand.grid(trait = 1:K, drug = 1:D)

  obs.stat <- list()
  U.complete <- list()
  V.complete <- list()
  for (d in 1:D) {
    U.temp <- V.temp <- list()
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
      if(test == "GEI") {
        U.bytrait[U.bytrait == 0] <- NA
        V.bytrait[V.bytrait == 0] <- NA
      }
      if(length(snp.list) == 1) {
        U.bytrait <- ifelse(is.na(diag(V.bytrait)), NA, U.bytrait)
      }
      obs.stat.temp[[k]] <- list(U = U.bytrait, V = V.bytrait)
      U.temp[[k]] <- U.bytrait
      V.temp[[k]] <- diag(V.bytrait)
    }
    obs.stat[[d]] <- obs.stat.temp
    U.complete[[d]] <- U.temp
    V.complete[[d]] <- V.temp
  }
  if(test == "GEI") {
    nonpolymorphic <- FALSE

    nonpoly.list <- list()
    analyze.snp <- list()
    for(k in 1:K) {
      ind <- NULL
      for(d in 1:D){
        ind <- rbind(ind, ifelse(is.na(obs.stat[[d]][[k]]$U), 1, 0))
      }
      nonpoly.list[[k]] <- which(apply(ind, 2, sum) > D-2)
      if(length(which(apply(ind, 2, sum)> D-2)) ==0){
        analyze.snp[[k]] <-  colnames(ind)
      }else{
        analyze.snp[[k]] <- colnames(ind)[-which(apply(ind, 2, sum)> D-2)]
      }
    }
    ## Remove non-polymorphic SNPs in at least two drug assignment groups
    # pre-process non-polymorphic SNPs that are polymorphic at least in two groups.
    if (any(is.na(unlist(U.complete)))) {

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
    if(all(is.na(unlist(V.complete))) | length(unlist(analyze.snp)) ==0){
      final.ret <- NA
      return(final.ret)
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
    Ucov <- list()
    for(iter in 1:nrow(trait.drug.ind)) {
      d1 <- trait.drug.ind[iter, 4]
      d2 <- trait.drug.ind[iter, 2]
      k1 <- trait.drug.ind[iter, 3]
      k2 <- trait.drug.ind[iter, 1]
      zeta.row <- k1 + (d1 - 1) * K
      zeta.col <- k2 + (d2 - 1) * K

      if (nonpolymorphic) {
        Ucov[[iter]] <- zeta.ret$zeta.GE[zeta.row, zeta.col] * diag(V.sqrt[[d1]][[k1]])  %*%
          t(diag(V.sqrt[[d2]][[k2]])) * update.R.C[[d1]][[which(trait.ind$row.ind == k1 & trait.ind$column.ind == k2)]]

      } else {
        Ucov[[iter]] <-  zeta.ret$zeta.GE[zeta.row, zeta.col] * diag(V.sqrt[[d1]][[k1]])  %*%
          t(diag(V.sqrt[[d2]][[k2]])) * R.C

      }

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
    beta.est <- NULL
    for(d in 1:D){
      for(k in 1:K) {
        tmp <- V.inv[[d]][[k]] %*% obs.stat1[[d]][[k]]$U[!is.na(obs.stat1[[d]][[k]]$U)]
        names(tmp) <- paste0(names(obs.stat1[[d]][[k]]$U), ":", k, ":",d)
        beta.est <- c(beta.est,tmp)
      }
    }
    trait.grid <- expand.grid(trait1 = 1:K, trait2 = 1:K)
    beta.cov <- list()
    for(d in 1:D) {
      beta.cov.tmp <- list()
      for(trait.pair in 1:nrow(trait.grid)) {
        k1 <- trait.grid[trait.pair, 1]
        k2 <- trait.grid[trait.pair, 2]
        if(k2 <= k1) {
          trait.drug.iter1 <- which(trait.drug.ind$row.trait == k1 & trait.drug.ind$column.trait == k2 &
                                      trait.drug.ind$row.drug == d & trait.drug.ind$column.drug == d)


          part1 <- V.inv[[d]][[k1]] %*%
            Ucov[[trait.drug.iter1]][!is.na(obs.stat1[[d]][[k1]]$U),
                                     !is.na(obs.stat1[[d]][[k2]]$U), drop = FALSE] %*%
            V.inv[[d]][[k2]]

          part1 <- part1[rownames(part1) %in% analyze.snp[[k1]],
                         colnames(part1) %in% analyze.snp[[k2]],drop = FALSE]
          part1.full <- matrix(NA, length(analyze.snp[[k1]]), length(analyze.snp[[k2]]))
          rownames(part1.full) <- analyze.snp[[k1]]
          colnames(part1.full) <- analyze.snp[[k2]]
          part1.full[match(rownames(part1), rownames(part1.full)),
                     match(colnames(part1), colnames(part1.full))] <- part1
        }else{
          part1.full <- t(beta.cov.tmp[[which(trait.grid$trait1 == k2 & trait.grid$trait2 == k1)]])
        }
        beta.cov.tmp[[trait.pair]] <- part1.full

      }
      names(beta.cov.tmp) <- apply(trait.grid, 1, function(x) paste0("Traits", x[1], "_", x[2]))
      beta.cov[[d]] <- beta.cov.tmp
    }

    ## prepare differentiated summary statistics for GEI test
    delta <- NULL
    for(d in setdiff(1:D, ref)) {
      for(k in 1:K) {
        delta.part1 <- delta.part2 <-  rep(NA, length(analyze.snp[[k]]))
        names(delta.part1) <- names(delta.part2) <- analyze.snp[[k]]

        delta.part1.tmp <- t(V.inv[[d]][[k]] %*% obs.stat1[[d]][[k]]$U[!is.na(obs.stat1[[d]][[k]]$U)])
        delta.part2.tmp <- t(V.inv[[ref]][[k]] %*% obs.stat1[[ref]][[k]]$U[!is.na(obs.stat1[[ref]][[k]]$U)])

        delta.part1[match(colnames(delta.part1.tmp), names(delta.part1))] <- delta.part1.tmp
        delta.part2[match(colnames(delta.part2.tmp), names(delta.part2))] <- delta.part2.tmp

        delta.diff <- delta.part1 - delta.part2
        names(delta.diff) <- paste0(analyze.snp[[k]], ":", k, ":", d)
        delta <- c(delta, delta.diff)
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
    rownames(Sigma.delta) <- colnames(Sigma.delta) <- names(delta)
    Sigma.delta <- Sigma.delta[!is.na(delta), !is.na(delta), drop = FALSE]

    final.ret <- list(beta.est = beta.est, beta.cov = beta.cov, delta = delta,
                      Sigma.delta = Sigma.delta, obs.stat = obs.stat,
                      nonpoly.list = nonpoly.list)
  }else if(test == "joint") {
    obs.stat1 <- obs.stat
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
      }
      U.cov <- NULL
      mlist <- list()
      for (col_id in 1:K) {
        diag.ind <- which(trait.ind$column.ind == col_id & trait.ind$row.ind == col_id)
        U.cov.col <- NULL
        index <- which(trait.ind$column.ind == col_id)
        row_id <- trait.ind[index, 1]
        for (iter in 1:length(index)) {
          U.cov.col <- rbind(U.cov.col, as.matrix(zeta.ret$zeta.joint[[d]][row_id[iter], col_id] * Ucov[[index[iter]]]))
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
      final.ret <- NA
      return(final.ret)
    }
    V.diag.all <- as.matrix(Matrix::bdiag(V.diag.bydrug))
    Sigma.inv <- V.diag.all %*% U.inv %*% V.diag.all #inverse of covariance matrix of beta

    U.all <- NULL
    for (DT.ind in 1:nrow(drug_trait)) {
      d <- drug_trait[DT.ind, 2]
      k <- drug_trait[DT.ind, 1]

      tmp <- obs.stat1[[d]][[k]]$U
      names(tmp) <- paste0(names(obs.stat1[[d]][[k]]$U), ":", k, ":",d)
      U.all <- c(U.all, tmp)
    }
    U.all <- as.matrix(U.all)
    beta.est <- MASS::ginv(V.diag.all) %*% U.all
    beta.cov <- MASS::ginv(Sigma.inv)
    rownames(beta.est) <- rownames(Sigma.inv) <- rownames(beta.cov) <-
      rownames(U.inv) <- colnames(U.inv) <-
      colnames(Sigma.inv) <- colnames(beta.cov) <- rownames(U.all)
    final.ret <- list(beta.est = beta.est, beta.cov = beta.cov,
                      U.all = U.all, Sigma.inv = Sigma.inv, V.diag.all = V.diag.all, U.inv = U.inv)
  }
  return(final.ret)
}

#' Calculate the summary statistics required by MTAGEI given the individual-level data
#'
#' This function allows you to calculate the summary statistics (beta, beta.cov) given the individual-level data.
#' @param dat a matrix containing traits and covariates data. Each row represent one subject.
#' @param geno.dat a matrix containing genotype data. Each row represent one subject.
#' @param cov.list a vector containing the covariates
#' @param SNP.list a vector containing the SNPs
#' @param trait.list a vector containing the traits
#' @param type type of the traits data.
#' @param ID indicator which column is unique subject ID
#' @param env which column represents the environmental variable
#' @param stats Use the score summary statistics or wald summary statistics. Default is "score".
#' @return summary statistics (beta, beta.cov) required by MTAGEI
#' @author Lan Luo
#' @export
#' @examples
#' \donttest{
#' data("rawdata")
#' names(rawdata)
#' attach(rawdata)
#'K <- 3
#'KA <- matrix(c(1, -0.03, 0.3, -0.03, 1, -0.535, 0.3, -0.535, 1), byrow = TRUE, nrow = 3)
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
#' ID = "ID", type = "continuous",
#' trait.list = paste0("Trait", 1:K),
#' SNP.list = colnames(geno.dat)[-1])
#' names(beta.sumstats.obj)
#' detach(rawdata)
#' }
Get_beta_cov_data <- function(dat, geno.dat,
                              cov.list, env, ID = "ID",
                              SNP.list, stats = "score",
                              trait.list, type){
  merge.dat <- merge(geno.dat, dat, by = ID)
  rownames(merge.dat) <- merge.dat$ID

  genotype <- merge.dat[, match(SNP.list, colnames(merge.dat)), drop = FALSE]
  geno.map <- data.frame(old = colnames(genotype), new = paste0("SNP", 1:ncol(genotype)))
  MAF <- colMeans(genotype)/2
  R.C <- cor(genotype)
  R.C[is.na(R.C)] <- 0
  colnames(genotype) <- geno.map$new
  simdata <- list(traits = merge.dat[, match(trait.list, colnames(merge.dat)), drop = FALSE],
                  covariates = merge.dat[, match(c(cov.list, env), colnames(merge.dat)), drop = FALSE])

  env.id <- grep(env, colnames(simdata$covariates))
  drug.list <- sort(unique(simdata$covariates[, env.id]))

  D <- length(drug.list) # number of different drug assignments
  K <- length(trait.list) # number of different drug response traits
  total.snp <- colnames(genotype)
  m.total <- length(total.snp) # number of rare variants in a gene
  L <- length(cov.list) # number of covariates

  MAC10.joint <- Reduce(union, sapply(1:D, function(x) colnames(genotype)[colSums(genotype[which(simdata$covariates[, env.id] == drug.list[x]), , drop = FALSE])>10], simplify = FALSE))
  MAC10.joint <- MAF[names(MAF) %in% geno.map[match(MAC10.joint, geno.map[, 2] ), 1]]
  highLD.joint <- try(caret::findCorrelation(R.C[colnames(R.C) %in% names(MAC10.joint), colnames(R.C) %in% names(MAC10.joint)],
                                             cutoff = 0.98),silent = T)
  if(length(highLD.joint)!=0 & !inherits(highLD.joint, "try-error")) {
    selSNP.joint<- names(MAC10.joint)[-highLD.joint]
  }else{
    selSNP.joint <- names(MAC10.joint)
  }

  MAC10.GEI <- Reduce(intersect, sapply(1:D, function(x) colnames(genotype)[colSums(genotype[which(simdata$covariates[, env.id] == drug.list[x]), , drop = FALSE])>10], simplify = FALSE))
  MAC10.GEI <- MAF[names(MAF) %in% geno.map[match(MAC10.GEI, geno.map[, 2] ), 1]]
  highLD.GEI <- try(caret::findCorrelation(R.C[colnames(R.C) %in% names(MAC10.GEI), colnames(R.C) %in% names(MAC10.GEI)],
                                           cutoff = 0.98),silent = T)
  if(length(highLD.GEI)!=0 & !inherits(highLD.GEI, "try-error")) {
    selSNP.GEI <- names(MAC10.GEI)[-highLD.GEI]
  }else{
    selSNP.GEI <- names(MAC10.GEI)
  }

  trait.grid <- expand.grid(trait1 = 1:K, trait2 = 1:K)
  beta.sumstats <- list()
  beta.sumstats.SPA <- list()
  beta.cov.list <- list()
  beta.est.list <- list()
  for(d in 1:D) {
    subject.ind <- which(simdata$covariates[, env.id] == drug.list[d])
    names(subject.ind) <- rownames(simdata$covariates)[subject.ind]

    X.alltraits.allSNP <- cbind(simdata$covariates[subject.ind, -env.id, drop = FALSE],
                                genotype[subject.ind, ,drop = FALSE])
    qr.X <- qr(X.alltraits.allSNP)
    X.alltraits <- X.alltraits.allSNP[, qr.X$pivot[1:qr.X$rank], drop = FALSE]#remove perfectly correlated SNPs

    dat.f <- as.formula(paste("Y", paste(colnames(X.alltraits), collapse = "+"), sep = "~") )
    dat.f.null <- as.formula(paste("Y", paste(colnames(X.alltraits)[-grep("SNP", colnames(X.alltraits))],
                                              collapse = "+"), sep = "~"))

    resid <- matrix(0, length(subject.ind), K)
    rownames(resid) <- names(subject.ind)
    beta.bytrait <- list()
    sigma.est <- NULL
    overlapping.samples <- list()
    covariates.list <- list()
    W.bytrait <- list()
    SPA <- ifelse(type == "binary", TRUE, FALSE)

    if(SPA) {
      V.SPA <- list()
      beta.sumstats.SPA.temp <- list()
      SPA.act <- list()
    }

    for(k in 1:K) {

      dat.full <- model.frame(dat.f, data.frame(Y = simdata$traits[subject.ind, k], X.alltraits))
      dat.full$id <- rownames(dat.full)

      if (any(is.na(dat.full))) {
        index <- na.omit(dat.full)$id
        dat.full <- dat.full[dat.full$id %in% index, ]
        m <- model.frame(dat.f, dat.full)
        mt <- attr(m, "terms")
        dat.full$response <- model.response(m, "numeric")
        mat <- as.data.frame(model.matrix(dat.f, m))
      }else {
        m <- model.frame(dat.f, dat.full)
        mt <- attr(m, "terms")
        dat.full$response <- model.response(m, "numeric")
        mat <- as.data.frame(model.matrix(dat.f, m))
      }
      qr.Xmat <- qr(mat)
      mat <- mat[, qr.Xmat$pivot[1:qr.Xmat$rank ], drop = FALSE] ## remove SNPs that linear dependent including the intercept term
      update.dat.f <- as.formula(paste("Y", paste(colnames(mat)[-1],
                                                  collapse = "+"), sep = "~") )
      update.dat.f.null <- as.formula(paste("Y", paste(colnames(mat)[-c(1, grep("SNP", colnames(mat)))],
                                                       collapse = "+"), sep = "~") )

      if(stats == "score") {
        ## make sure the genotype matrix is full rank.
        geno.dat <- as.matrix(mat[, grep("SNP", colnames(mat)), drop = FALSE])
        cov.dat <- as.matrix((mat[, -grep("SNP", colnames(mat)), drop = FALSE]))
      }

      if(type == "continuous") {

        if(stats == "score") {
          mod1 <- lm(update.dat.f.null, data = dat.full)

          ## score statistics
          s2.null <- summary(mod1)$sigma^2
          resid.null <- mod1$residuals

          U <- as.vector(crossprod(geno.dat, resid.null)) / s2.null
          names(U) <- geno.map[match(colnames(geno.dat), geno.map[, 2]), 1]
          W <- crossprod(geno.dat, geno.dat) - crossprod(geno.dat, cov.dat) %*%
            solve(crossprod(cov.dat, cov.dat)) %*% crossprod(cov.dat, geno.dat)
          V <- W / s2.null
          if(m.total > 2) {
            colnames(V) <- rownames(V) <- geno.map[match(colnames(geno.dat), geno.map[, 2]), 1]
          }else{
            names(V) <- geno.map[match(colnames(geno.dat), geno.map[, 2]), 1]
          }
          beta.est <- as.vector(U %*% MASS::ginv(V))
          names(beta.est) <- geno.map[match(colnames(geno.dat), geno.map[, 2]), 1]
          beta.temp <- beta.est
          covariate.list.temp <- colnames(mat)
          residual.temp <- mod1$residuals

        }else if(stats == "wald") {
          mod1 <- lm(update.dat.f, data = dat.full)
          beta.est <- coef(mod1)
          beta.est <- beta.est[!is.na(beta.est)]
          beta.temp <-  beta.est[grep("SNP", names(beta.est))]
          covariate.list.temp <- names(beta.est)
          names(beta.temp) <- geno.map[match(names(beta.temp), geno.map[, 2]), 1]
          residual.temp <- mod1$residuals

        }
        sigma.est <- c(sigma.est, summary(mod1)$sigma^2)

      }else if(type == "binary") {
        mod1 <- glm(update.dat.f.null, data = dat.full, family = binomial(logit))
        if(stats == "score") {
          b1 <- mod1$fitted.values #equivalent to mod1$predict = exp(mod1$linear.predictors)/(1 + exp(mod1$linear.predictors))

          b2 <- b1 * (1 - b1)
          b2.prod <- b2 %*% matrix(1, 1, ncol(cov.dat)) * cov.dat
          #### b2.prod is equivalent to t(cov.dat) %*% diag(as.vector(b2), length(b2)) %*% cov.dat,since Sigma.inv is diagonal matrix

          U <- crossprod(geno.dat, (dat.full$Y - b1))
          U <- as.vector(U)
          names(U) <- geno.map[match(colnames(geno.dat), geno.map[, 2]), 1]
          Sigma.inv <- diag(as.vector(b2), length(b2))
          Sigma.cov <- b2 * cov.dat  #crossprod(Sigma.inv, cov.dat)

          V.tmp <- Sigma.inv - Sigma.cov %*% tcrossprod(solve(crossprod(cov.dat, Sigma.cov)), Sigma.cov)
          V <- crossprod(geno.dat, V.tmp) %*% geno.dat
          if(m.total > 2) {
            colnames(V) <- rownames(V) <- geno.map[match(colnames(geno.dat), geno.map[, 2]), 1]
          }else{
            names(V) <- geno.map[match(colnames(geno.dat), geno.map[, 2]), 1]
          }
          beta.est <- as.vector(U %*% MASS::ginv(V))
          names(beta.est) <- geno.map[match(colnames(geno.dat), geno.map[, 2]), 1]
          beta.temp <- beta.est
          covariate.list.temp <- colnames(mat)
          residual.temp <- (dat.full$Y- b1) / sqrt(b2)
          W.bytrait[[k]] <-  diag(b2)

          if(SPA) {
            ## SPA adjustment
            SPA.ret <- SPAtest::ScoreTest_SPA_wMeta(genos = t(geno.dat), pheno = dat.full$Y, cov = cov.dat,
                                                    method = "fastSPA", output = "metaZ", minmac = 10)

            SPA.index <- (SPA.ret$p.value != SPA.ret$p.value.NA)
            names(SPA.index) <- geno.map[match(colnames(geno.dat), geno.map[, 2]), 1]
            SPA.act[[k]] <- SPA.index

            Z <- sign(SPA.ret$p.value) * qnorm(abs(SPA.ret$p.value)/2, lower.tail = FALSE)
            adj.SNP <- colnames(geno.dat)[!is.na(Z)]
            names(Z) <- geno.map[match(colnames(geno.dat), geno.map[, 2]), 1]
            beta.SPA.temp <- Z^2/U
            beta.sumstats.SPA.temp[[k]] <- beta.SPA.temp
            V.SPA.diag <- (U/Z)^2
            V.SPA[[k]] <- V.SPA.diag
          }
        }
      }

      resid[rownames(resid) %in% rownames(mat), k] <- residual.temp
      beta.bytrait[[k]] <- beta.temp
      covariates.list[[k]] <- covariate.list.temp
      overlapping.samples[[k]] <- rownames(mat)
    }
    names(beta.bytrait) <-  paste0("Trait", 1:K)
    len <- count.overlapping(overlapping.samples) # Pan
    len[len < 0] <- 0

    ave.rr <- t(resid) %*% resid/(len) ## updated 1/N sum rr^T
    ave.rr[is.nan(ave.rr)] <- 0
    betacov.bytrait <- list()
    betacov.bytrait.SPA <- list()

    for(trait.pair in 1:nrow(trait.grid)) {
      k1 <- trait.grid[trait.pair, 1]
      k2 <- trait.grid[trait.pair, 2]
      if(k2 <= k1) {
        X1 <- as.matrix(cbind(1, X.alltraits[rownames(X.alltraits) %in% overlapping.samples[[k1]],
                                             colnames(X.alltraits) %in% covariates.list[[k1]], drop = FALSE]))
        X2 <- as.matrix(cbind(1, X.alltraits[rownames(X.alltraits) %in% overlapping.samples[[k2]],
                                             colnames(X.alltraits) %in% covariates.list[[k2]], drop = FALSE]))
        overlap.subset <- intersect(rownames(X1), rownames(X2))
        if(type == "continuous") {
          step11 <- crossprod(X1, X1) #t(X1) %*% X1
          step22 <- crossprod(X2, X2) #t(X2) %*% X2
          step12 <- crossprod((X1[rownames(X1) %in% overlap.subset, , drop = FALSE]),
                              X2[rownames(X2) %in% overlap.subset, , drop = FALSE])
        }else if (type == "binary") {
          step11 <- crossprod(X1 * diag(W.bytrait[[k1]]), X1) # t(X1) %*% W.bytrait[[k1]] %*% X1
          step22 <- crossprod(X2 * diag(W.bytrait[[k2]]), X2) # t(X2) %*% W.bytrait[[k2]] %*% X2

          W.sqrt1 <- sqrt(diag(W.bytrait[[k1]]))
          W.sqrt2 <- sqrt(diag(W.bytrait[[k2]]))
          step12 <- crossprod(X1[rownames(X1) %in% overlap.subset, , drop = FALSE] * W.sqrt1, W.sqrt2 *
                                X2[rownames(X2) %in% overlap.subset, , drop = FALSE])

        }
        betacov.tmp <- solve(step11) %*% (ave.rr[k1, k2] * step12) %*% solve(step22)

        betacov.tmp <- betacov.tmp[grep("SNP", rownames(betacov.tmp)),
                                   grep("SNP", colnames(betacov.tmp)), drop = FALSE]
        rownames(betacov.tmp) <- geno.map[match(rownames(betacov.tmp), geno.map[, 2] ), 1]
        colnames(betacov.tmp) <- geno.map[match(colnames(betacov.tmp), geno.map[, 2] ), 1]


      }else{
        betacov.tmp <- t(betacov.bytrait[[which(trait.grid$trait1 == k2 & trait.grid$trait2 == k1)]])

      }

      betacov.bytrait[[trait.pair]] <- betacov.tmp

    }
    names(betacov.bytrait) <- apply(trait.grid, 1, function(x) paste0("Traits", x[1], "_", x[2]))
    beta.sumstats[[d]] <- list(beta.est = beta.bytrait, beta.cov = betacov.bytrait)
    if(SPA) {
      beta.sumstats.SPA[[d]] <- list(beta.est = beta.sumstats.SPA.temp,
                                     # beta.cov = betacov.bytrait.SPA,
                                     V.SPA = V.SPA,
                                     SPA.act = SPA.act)
    }

    beta.est.list[[d]] <- beta.sumstats[[d]]$beta.est
    beta.cov.list[[d]] <- beta.sumstats[[d]]$beta.cov
  }


  final.ret <- list(beta.est = beta.est.list,
                    beta.cov = beta.cov.list,
                    beta.sumstats = beta.sumstats,
                    beta.sumstats.SPA = beta.sumstats.SPA,
                    selSNP.joint = selSNP.joint,
                    selSNP.GEI = selSNP.GEI)
  # names(beta.sumstats) <- paste0("Drug", drug.list)
  # return(beta.sumstats)
}

#' Calculate Covariances of Z-scores between Traits from Overlapping Samples
#'
#' This function allows you to estimate the matrix \eqn{\zeta} to adjust for the potential sample overlap in the data set. Here we applied LD pruning (\eqn{r^2 < 0.1} in 500kb region) on 1000 Genome genotype dataset (hg19) as a list of reference independent SNPs. The SNP ID is chr:pos.
#' @param common.sumstats a numeric list of score summary statistics for all independent common variants that are not associated with the traits
#' @param trait.list a vector containing the traits
#' @return A \eqn{K \times K} matrix \eqn{\zeta}, where \eqn{K} is the number of traits.
#' @author Lan Luo
#' @export
#' @examples
#' \donttest{data(zeta.example)
#' data("sumstats.dat")
#' names(sumstats.dat)
#' attach(sumstats.dat)
#' str(common.sumstats)
#' K <- 3
#' zeta.ret <- Get_zeta(common.sumstats = common.sumstats, trait.list = paste("Trait", 1:K))
#' zeta.ret
#' detach(sumstats.dat)
#' }
Get_zeta <- function(common.sumstats, trait.list) {
  drug.list <- names(common.sumstats)
  D <- length(drug.list) # number of different drug assignments
  K <- length(trait.list) # number of different drug response traits
  m.common <- length(common.sumstats[[1]][[1]]$U)

  message(paste0("Begin the calculation of summary statistics for ", m.common, " common independent null SNPs"))
  zeta.list <- list()
  for(d in 1:D) {
    common.zscore <- matrix(0, m.common, K)
    for(k in 1:K) {
      common.zscore[, k] <- common.sumstats[[d]][[k]]$U/ sqrt(common.sumstats[[d]][[k]]$V)
    }
    pval <- 2 * pnorm(-abs(common.zscore))
    ind <- apply(pval, 1, function(x) ifelse(any(x < 0.05), 1, 0))
    zeta.temp <-cor(common.zscore[ind == 0, ])
    rownames(zeta.temp) <- colnames(zeta.temp) <- paste0("Trait", 1:K)
    zeta.list[[d]] <- zeta.temp
  }
  names(zeta.list) <- paste0("grp", 1:D)

  common.sumstats.combined <- list()
  for(k in 1:K) {
    U.temp1 <- common.sumstats[[1]][[k]]$U
    V.temp1 <- common.sumstats[[1]][[k]]$V
    for(d in 2:D) {
      U.temp1 <- U.temp1 + common.sumstats[[d]][[k]]$U
      V.temp1 <- V.temp1 + common.sumstats[[d]][[k]]$V
    }
    common.sumstats.combined[[k]] <- list(U = U.temp1, V = V.temp1)
  }

  common.zscore <- matrix(0, m.common, K)
  for(k in 1:K) {
    common.zscore[, k] <- common.sumstats.combined[[k]]$U/ sqrt(common.sumstats.combined[[k]]$V)
  }
  pval <- 2 * pnorm(-abs(common.zscore))
  ind <- apply(pval, 1, function(x) ifelse(any(x < 0.05), 1, 0))
  zeta <- cor(common.zscore[ind == 0, ])
  rownames(zeta) <- colnames(zeta) <- paste0("Trait", 1:K)

  zeta.GE <- matrix(0, D*K, D*K)
  for(d in 1:D) {
    ind <- 1:K + (d - 1) * K
    zeta.GE[ind, ind] <- zeta.list[[d]]
  }
  rownames(zeta.GE) <- colnames(zeta.GE) <- as.vector(sapply(1:D, function(x) paste0(paste0("grp", x),  paste0("Trait", 1:K))))
  zeta.ret <- list(zeta.main = zeta,
                   zeta.joint = zeta.list,
                   zeta.GE = zeta.GE)
  return(zeta.ret)
}
