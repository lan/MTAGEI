# non-exported functions from SKAT package
# https://github.com/leeshawn/SKAT
# date: 01/23/2020

Get_Lambda<-function (K)
{
  out.s <- eigen(K, symmetric = TRUE, only.values = TRUE)
  lambda1 <- out.s$values
  IDX1 <- which(lambda1 >= 0)
  IDX2 <- which(lambda1 > mean(lambda1[IDX1])/1e+05)
  if (length(IDX2) == 0) {
    stop("No Eigenvalue is bigger than 0!!")
  }
  lambda <- lambda1[IDX2]
  return(lambda)
}
Get_Liu_PVal.MOD.Lambda <- function (Q.all, lambda, log.p = FALSE)
{
  param <- Get_Liu_Params_Mod_Lambda(lambda)
  Q.Norm <- (Q.all - param$muQ)/param$sigmaQ
  Q.Norm1 <- Q.Norm * param$sigmaX + param$muX
  p.value <- pchisq(Q.Norm1, df = param$l, ncp = param$d,
                    lower.tail = FALSE, log.p = log.p)
  return(p.value)
}

Get_Liu_Params_Mod_Lambda <- function (lambda)
{
  c1 <- rep(0, 4)
  for (i in 1:4) {
    c1[i] <- sum(lambda^i)
  }
  muQ <- c1[1]
  sigmaQ <- sqrt(2 * c1[2])
  s1 = c1[3]/c1[2]^(3/2)
  s2 = c1[4]/c1[2]^2
  beta1 <- sqrt(8) * s1
  beta2 <- 12 * s2
  type1 <- 0
  if (s1^2 > s2) {
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 * a^3 - a^2
    l = a^2 - 2 * d
  }
  else {
    type1 <- 1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <- l + d
  sigmaX <- sqrt(2) * a
  re <- list(l = l, d = d, muQ = muQ, muX = muX, sigmaQ = sigmaQ,
             sigmaX = sigmaX)
  return(re)
}

Get_Liu_PVal.MOD.Lambda.Zero <- function (Q, muQ, muX, sigmaQ, sigmaX, l, d)
{
  Q.Norm <- (Q - muQ)/sigmaQ
  Q.Norm1 <- Q.Norm * sigmaX + muX
  temp <- c(0.05, 10^-10, 10^-20, 10^-30, 10^-40, 10^-50,
            10^-60, 10^-70, 10^-80, 10^-90, 10^-100)
  out <- qchisq(temp, df = l, ncp = d, lower.tail = FALSE)
  IDX <- max(which(out < Q.Norm1))
  pval.msg <- sprintf("Pvalue < %e", temp[IDX])
  return(pval.msg)
}

Get_Liu_Params_Mod <- function (c1)
{
  muQ <- c1[1]
  sigmaQ <- sqrt(2 * c1[2])
  s1 = c1[3]/c1[2]^(3/2)
  s2 = c1[4]/c1[2]^2
  beta1 <- sqrt(8) * s1
  beta2 <- 12 * s2
  type1 <- 0
  if (s1^2 > s2) {
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 * a^3 - a^2
    l = a^2 - 2 * d
  }
  else {
    type1 <- 1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <- l + d
  sigmaX <- sqrt(2) * a
  re <- list(l = l, d = d, muQ = muQ, muX = muX, sigmaQ = sigmaQ,
             sigmaX = sigmaX)
  return(re)
}

SKAT_Optimal_Integrate_Func_Davies <- function (x, pmin.q, param.m, r.all)
{
  n.r <- length(r.all)
  n.x <- length(x)
  temp1 <- param.m$tau %x% t(x)
  temp <- (pmin.q - temp1)/(1 - r.all)
  temp.min <- apply(temp, 2, min)
  re <- rep(0, length(x))
  for (i in 1:length(x)) {
    min1 <- temp.min[i]
    if (min1 > sum(param.m$lambda) * 10^4) {
      temp <- 0
    }
    else {
      min1.temp <- min1 - param.m$MuQ
      sd1 <- sqrt(param.m$VarQ - param.m$VarRemain)/sqrt(param.m$VarQ)
      min1.st <- min1.temp * sd1 + param.m$MuQ
      dav.re <- CompQuadForm::davies(min1.st, param.m$lambda, acc = 10^(-6))
      temp <- dav.re$Qq
      if (dav.re$ifault != 0) {
        stop("dav.re$ifault is not 0")
      }
    }
    if (temp > 1) {
      temp = 1
    }
    re[i] <- (1 - temp) * dchisq(x[i], df = 1)
  }
  return(re)
}
SKAT_Optimal_PValue_Liu <- function (pmin.q, param.m, r.all, pmin = NULL)
{
  re <- integrate(SKAT_Optimal_Integrate_Func_Liu, lower = 0,
                  upper = 40, subdivisions = 2000, pmin.q = pmin.q, param.m = param.m,
                  r.all = r.all, abs.tol = 10^-25)
  pvalue <- 1 - re[[1]]
  if (!is.null(pmin)) {
    if (pmin * length(r.all) < pvalue) {
      pvalue = pmin * length(r.all)
    }
  }
  return(pvalue)
}
SKAT_Optimal_Integrate_Func_Liu <- function (x, pmin.q, param.m, r.all)
{
  n.r <- length(r.all)
  n.x <- length(x)
  temp1 <- param.m$tau %x% t(x)
  temp <- (pmin.q - temp1)/(1 - r.all)
  temp.min <- apply(temp, 2, min)
  temp.q <- (temp.min - param.m$MuQ)/sqrt(param.m$VarQ) *
    sqrt(2 * param.m$Df) + param.m$Df
  re <- pchisq(temp.q, df = param.m$Df) * dchisq(x, df = 1)
  return(re)
}

Get_Liu_Params <- function (c1)
{
  muQ <- c1[1]
  sigmaQ <- sqrt(2 * c1[2])
  s1 = c1[3]/c1[2]^(3/2)
  s2 = c1[4]/c1[2]^2
  beta1 <- sqrt(8) * s1
  beta2 <- 12 * s2
  type1 <- 0
  if (s1^2 > s2) {
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 * a^3 - a^2
    l = a^2 - 2 * d
  }
  else {
    type1 <- 1
    a = 1/s1
    d = 0
    l = 1/s1^2
  }
  muX <- l + d
  sigmaX <- sqrt(2) * a
  re <- list(l = l, d = d, muQ = muQ, muX = muX, sigmaQ = sigmaQ,
             sigmaX = sigmaX)
  return(re)
}

Get_Liu_PVal <- function (Q, W, Q.resampling = NULL)
{
  Q.all <- c(Q, Q.resampling)
  A1 <- W/2
  A2 <- A1 %*% A1
  c1 <- rep(0, 4)
  c1[1] <- sum(diag(A1))
  c1[2] <- sum(diag(A2))
  c1[3] <- sum(A1 * t(A2))
  c1[4] <- sum(A2 * t(A2))
  param <- Get_Liu_Params(c1)
  Q.Norm <- (Q.all - param$muQ)/param$sigmaQ
  Q.Norm1 <- Q.Norm * param$sigmaX + param$muX
  p.value <- pchisq(Q.Norm1, df = param$l, ncp = param$d,
                    lower.tail = FALSE)
  p.value.resampling = NULL
  if (length(Q.resampling) > 0) {
    p.value.resampling <- p.value[-1]
  }
  re <- list(p.value = p.value[1], param = param, p.value.resampling = p.value.resampling)
  return(re)
}

Get_PValue.Lambda <- function (lambda, Q)
{
  n1 <- length(Q)
  p.val <- rep(0, n1)
  p.val.liu <- rep(0, n1)
  is_converge <- rep(0, n1)
  p.val.liu <- Get_Liu_PVal.MOD.Lambda(Q, lambda)
  for (i in 1:n1) {
    out <- CompQuadForm::davies(Q[i], lambda, acc = 10^(-6))
    p.val[i] <- out$Qq
    is_converge[i] <- 1
    if (length(lambda) == 1) {
      p.val[i] <- p.val.liu[i]
    }
    else if (out$ifault != 0) {
      is_converge[i] <- 0
    }
    if (p.val[i] > 1 || p.val[i] <= 0) {
      is_converge[i] <- 0
      p.val[i] <- p.val.liu[i]
    }
  }
  p.val.msg = NULL
  p.val.log = NULL
  if (p.val[1] == 0) {
    param <- Get_Liu_Params_Mod_Lambda(lambda)
    p.val.msg <- Get_Liu_PVal.MOD.Lambda.Zero(Q[1], param$muQ,
                                              param$muX, param$sigmaQ, param$sigmaX, param$l,
                                              param$d)
    p.val.log <- Get_Liu_PVal.MOD.Lambda(Q[1], lambda, log.p = TRUE)[1]
  }
  return(list(p.value = p.val, p.val.liu = p.val.liu, is_converge = is_converge,
              p.val.log = p.val.log, pval.zero.msg = p.val.msg))
}
