#' Approach 2:
#' IPW Estimation and Hypothesis Testing of Strain-Specific
#' Vaccine Efficacy Adjusted for Covariate Effects with Missing Causes
#'
#' Fit a stratified cause-specific proportional hazards regression model
#' postulates that the conditional cause-specific hazard function
#' for cause j for an individual in the k-th stratum with the covariate value
#' \eqn{z} equals
#' \deqn{ \lambda_{kj}(t|z)=\lambda_{0kj}(t)\exp(\beta_j^{T}z) }
#' for j=1,...,J and k=1, ..., K, where \eqn{\lambda_{0kj}(\cdot)}
#' is an unspecified cause-specific baseline hazard function
#' for the k-th stratum and K is the number of strata.
#'
#' Approach 2: Failure endpoint = COVID.
#' Cause/mark V=0: vaccine-matched genotype AND viral load above minimum threshold;
#' Cause/mark V=1: vaccine-mismatched genotype AND viral load above minimum threshold;
#' Cause/mark V=2: viral load below minimum threshold.
#' In Approach 2, we only care about V=0 and V=1;
#' V=2 is a competing risk for both V=0 and V=1.
#'
#' Let R be the indicator of whether all possible data are observed for
#' a subject; R = 1 if either \eqn{\delta=0} (right-censored) or if
#' \eqn{\delta=1} (a failure is observed) and the cause is not missing;
#' and R = 0 otherwise.
#' Let VL be the viral load and \eqn{h_0} the minimum threshold.
#' Missing at random (MAR) assumption is valid using Approach 2:
#' \eqn{P(R=1|\delta=1,T,Z,VL,h_0,V)=P(R=1|\delta=1,T,Z,VL,h_0)}
#' \eqn{=I(VL<h_0)+P(R=1|\delta=1,T,Z,VL,VL \ge h_0)I(VL\ge h_0)}
#'
#'
#' @param cmprskPHformula a formula object with the response on the left of a '~'
#' operator, and the independent terms on the right as covariates in the
#' cause-specific proportional hazards regression model. The response
#' must be in the format of cbind(time, delta, cause), where time is the minimum
#' of the failure time and the censoring time, delta specifies the censoring status
#' (1: failure is observed; 0: right-censored), cause is the cause of failure.
#' If cause is missing, cause=NA.
#' @param trtpos the position of the treatment group indicator in the RHS of
#' the formula object cmprskPHformula. The default value is 1.
#' @param strata the column name of the strata variable in the data.
#' @param causelevels types of causes.
#' @param missmodel indicates whether a logistic regression model is applied to
#' predict the probability of being not missing for each type of causes under the
#' missing at random (MAR) assumption. If not, a logical vector is needed. TRUE:
#' fit a logistic regression model; FALSE: that cause is not missing, i.e. the
#' probability of being not missing for that cause is 1. For example, we have
#' three types of causes: c("cause A", "cause B", "cause C"). missmodel=c(T,T,F)
#' means that the probability of being not missing for cause C is 1, and the
#' logistic regression model is fitted for cause A and cause B.
#' The default value is T.
#' @param missformula a formula object started with a '~' operator, and the independent
#' terms on the right of '~' are variables used for predicting the probability
#' of being not missing under a logistic regression model.
#' @param data a data.frame with the variables.
#' @param VEnull the assumed VE value in the null hypothesis. The default value
#' is 0.3.
#' @param maxit Maximum number of iterations to attempt for convergence. The
#' default is 15.
#'
#' @return returns an object of type 'markPH.ipw2'. With the following arguments:
#' \item{causes}{types of causes}
#' \item{misscauses}{types of causes which are subject to missingness}
#' \item{coef.cc}{estimates of covariate coefficients using a complete-case analysis}
#' \item{se.cc}{estimates of standard errors of estimators for covariate coefficients
#' using a complete-case analysis}
#' \item{coef}{estimates of covariate coefficients using IPW method}
#' \item{se}{estimates of standard errors of estimators for covariate coefficients
#' using IPW method}
#' \item{coef.VE}{estimates of the strain-specific vaccine efficacy to reduce
#' susceptibility to strain j using IPW method, j=1,...,J}
#' \item{se.VE}{estimates of standard errors of estimators for the strain-specific
#' vaccine efficacy to reduce susceptibility to strain j using IPW method, j=1,...,J}
#' \item{coef.VD}{estimates of VD using IPW method}
#' \item{se.VD}{estimates of standard errors of estimators for VD using IPW method}
#' \item{U1}{test statistic for testing against the alternative hypothesis
#' HA1: VE(j)>=VEnull with strict inequality for some j in \code{misscauses}}
#' \item{U2}{test statistic for testing against the alternative hypothesis
#' HA2: VE(j) is not equal to VEnull for some j in \code{misscauses}}
#' \item{U1j}{test statistic for testing against the alternative hypothesis
#' HAj1: VE(j)>VEnull}
#' \item{U2j}{test statistic for testing against the alternative hypothesis
#' HAj2: VE(j) is not equal to VEnull}
#' \item{T1}{test statistic for testing against the alternative hypothesis
#' HB1: VE(\code{misscauses}[1]) \eqn{\ge \dots \ge} VE(\code{misscauses}[I]) with at least one strict inequality}
#' \item{T2}{test statistic for testing against the alternative hypothesis
#' HB2: VE(i) is not equal to VE(j) for at least one pair of i and j in \code{misscauses}}
#' \item{pval.A1}{p value for testing against the alternative hypothesis
#' HA1: VE(j)>=VEnull with strict inequality for some j in \code{misscauses}}
#' \item{pval.A2}{p value for testing against the alternative hypothesis
#' HA2: VE(j) is not equal to VEnull for some j in \code{misscauses}}
#' \item{pval.A1j}{p value for testing against the alternative hypothesis
#' HAj1: VE(j)>VEnull}
#' \item{pval.A2j}{p value for testing against the alternative hypothesis
#' HAj2: VE(j) is not equal to VEnull}
#' \item{pval.B1}{p value for testing against the alternative hypothesis
#' HB1: VE(\code{misscauses}[1]) \eqn{\ge \dots \ge} VE(\code{misscauses}[I]) with at least one strict inequality}
#' \item{pval.B2}{p value for testing against the alternative hypothesis
#' HB2: VE(i) is not equal to VE(j) for at least one pair of i and j in \code{misscauses}}
#' \item{tgrid}{specifies the time grids for estimating \eqn{\lambda_{0kj}(\cdot)}
#' , the unspecified cause-specific baseline hazard function
#' for cause j and the k-th stratum.}
#' \item{lambda0}{estimates of \eqn{\lambda_{0kj}(\cdot)} using IPW method}
#' \item{covariates}{covariates in the cause-specific proportional hazards model.}
#' \item{trtpos}{the position of the treatment group indicator in the RHS of
#' the formula object cmprskPHformula.}
#' \item{VEnull}{the assumed VE value in the null hypothesis.}
#' \item{diffsigma}{estimates of standard deviation of \eqn{\hat\alpha_j-\hat\alpha_{j-1}}, j=2,...,J.}
#' \item{cov.alpha}{estimated covariance matrix of \eqn{(\hat\alpha_1,\dots,\hat\alpha_J)}.}
#'
#' @author Fei Heng
#'
#' @references (2021+)
#'
#' @examples
#'
#' ## Example 1: Simulated competing risks data with 2 causes subject to missingness
#' data(sim2cs)
#'
#' threshold <- 0.5
#' sim2cs$cause[(sim2cs$delta==1)&(sim2cs$A<threshold)] <- 3
#'
#' res.ipw2 <- markPH.ipw2(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
#'                         trtpos=1,
#'                         strata="strata",
#'                         causelevels=c(1,2,3),
#'                         missmodel=c(TRUE,TRUE,FALSE),
#'                         missformula=~z1+A,
#'                         data=sim2cs,
#'                         VEnull=0.3,
#'                         maxit=15)
#' res.ipw2 # print the result
#'
#' ## Example 2: Simulated competing risks data with 3 causes subject to missingness
#' data(sim3cs)
#'
#' threshold <- 0.5
#' sim3cs$cause[(sim3cs$delta==1)&(sim3cs$A<threshold)] <- 4
#'
#' res.ipw2 <- markPH.ipw2(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
#'                         trtpos=1,
#'                         strata="strata",
#'                         causelevels=c(1,2,3,4),
#'                         missmodel=c(TRUE,TRUE,TRUE,FALSE),
#'                         missformula=~z1+A,
#'                         data=sim3cs,
#'                         VEnull=0.3,
#'                         maxit=15)
#' res.ipw2 # print the result
#'
#' @import MASS nnet
#'
#' @export



markPH.ipw2 <- function(cmprskPHformula,
                       trtpos=1,
                       strata,
                       causelevels=NULL,
                       missmodel=T,
                       missformula,
                       data=parent.frame(),
                       VEnull=0.3,
                       maxit=15){
  # next version:
  # ...

  if (all(missmodel)) {
    stop("Please use markPH.ipw().")
  } else if (all(!missmodel)) {
    stop("Please use markPH().")
  }

  data[data=="NaN"] <- NA

  a <- model.frame(cmprskPHformula, data=data, na.action = na.pass)
  time <- model.response(a)[,1]
  delta <- model.response(a)[,2] # 0: censored, 1: observed failure
  if (is.null(causelevels)){
    cause.fa <- factor(model.response(a)[,3])
  } else {
    cause.fa <- factor(model.response(a)[,3], labels=causelevels)
  }
  causelevels <- levels(cause.fa)
  cause <- as.numeric(cause.fa) # cause = 1,2,...,J
  cause[is.na(cause.fa)] <- NA


  covar2 <- as.matrix(model.matrix(a, data=a, na.action = na.pass)[,-1])
  #A <- data[,Aux]
  strata.fa <- factor(data[,strata])
  strata.num <- as.numeric(strata.fa)

  R <- (delta==0)|(!is.na(cause))
  data$R <- R # cen.code=0

  nstrt <- length(levels(strata.fa))# number of strata
  ncs <- length(unique(cause[!is.na(cause)]))# number of causes
  ncov <- ncol(covar2)# number of covariates
  nsamp <- length(time)


  # set grid points for lambda0
  time0 <- min(time)
  tau <- max(time)
  ngrid <- 50 # number of grids (hard coding)
  tgrid <- seq(time0,tau,length.out=ngrid)
  b <- (tau-time0)/5 # bandwidth for lambda_0(t) (hard coding)




  missformula <- formula(paste(c("R",missformula),collapse = " "))
  ## predict P(R=1) in jj-th stratum
  # logistic regression
  wipw <- R
  a <- model.frame(missformula, data=data, na.action = na.pass)
  covar.miss <- as.matrix(model.matrix(a, data=a, na.action = na.pass)[,-1])

  npsi <- ncol(covar.miss)+1
  dr <- matrix(0,nsamp,npsi)
  Ipsi <- array(0, c(npsi,npsi,nstrt))
  Spsi <- matrix(0,nsamp,npsi)
  for (jj in 1:nstrt){

    temp <- which((delta==1)&(strata.num==jj)&!(cause.fa %in% causelevels[!missmodel]))
    miss.res <- glm(missformula, data=data, family='binomial',subset=temp)
    rhat <- predict(miss.res, as.data.frame(covar.miss[temp,]), type="response")

    Rjj <- R[temp]
    njj <- length(Rjj)
    drjj <- matrix(0,njj,npsi)
    Spsijj <- matrix(0,njj,npsi)

    W <- cbind(rep(1,njj), covar.miss[temp,])

    for (ipsi in 1:npsi){
      Spsijj[,ipsi] <- (Rjj-rhat)*W[,ipsi]
      drjj[,ipsi] <- rhat*(1-rhat)*W[,ipsi]
    }
    Ipsi[,,jj] <- crossprod(W,drjj)/njj

    dr[temp,] <- drjj
    Spsi[temp,] <- Spsijj
    wipw[temp] <- Rjj/rhat
  }

  # max(wipw[cause.fa %in% causelevels[missmodel]])


  ## estimation for each cause
  # initialization
  sbeta_c <- matrix(0, ncov, ncs)
  sstd_c <- matrix(0, ncov, ncs)

  sbeta_ic <- matrix(0, ncov, ncs)
  sstd_ic <- matrix(0, ncov, ncs)

  sVE_ic <- rep(0,ncs)
  sVEstd_ic <- rep(0,ncs)
  sVD_ic <- rep(0,ncs-1)
  sVDstd_ic <- rep(0,ncs-1)


  for (ics in 1:ncs){
    # deltacs: indicator of non-censoring for each cause
    deltacs <- cause == ics# be carefule: NaN ->F
    deltacs[is.na(cause)] <- F

    # complete estimation
    beta0 <- rep(0,ncov)
    res.cc <- estf(time,covar2,deltacs,beta0,strata.num,maxit,subset=R)

    # ipw estimation
    beta0 <- rep(0,ncov)
    res.ipw <- esti(time,covar2,deltacs,beta0,wipw,strata.num,maxit,nstrt,dr,Ipsi,Spsi)


    # complete result
    sbeta_c[,ics] <- res.cc$est
    sstd_c[,ics] <- res.cc$se
    # ipw-c result
    sbeta_ic[,ics] <- res.ipw$est
    sstd_ic[,ics] <- res.ipw$se
  }

  ## calculate the cumulative baseline function
  # ipw
  # sLamtA0_c <- array(0, c(ngrid,nstrt,ncs))
  slamtA0_c <- array(0, c(ngrid,nstrt,ncs))
  # slamtA0_c_bd <- array(0, c(ngrid,nstrt,ncs)) boundary correction
  for (ics in 1:ncs){
    deltacs <- cause == ics# be careful: NaN ->F
    deltacs[is.na(cause)] <- F


    ebz <- exp(covar2%*%sbeta_ic[,ics])
    wt <- wipw*deltacs

    Sf_ic <- rep(0, nsamp)
    for (i in 1:nsamp){
      if (wt[i]>0){
        Sf_ic[i] <- sum(ebz*wipw*(time >= time[i])*(strata.num==strata.num[i]))
      }
    }

    for (jj in 1:nstrt){
      for (it in 1:ngrid){
        tvalue <- tgrid[it]

        # temp=which((wt!=0)&(Sf_ic!=0)&(time<tvalue)&(strata.num==jj))
        # sLamtA0_c[it,jj,ics]  <- sum(wt[temp]/Sf_ic[temp])

        temp=which((Sf_ic!=0)&(strata.num==jj))
        slamtA0_c[it,jj,ics] <- sum(epanker(tvalue, time[temp], b)*wt[temp]/Sf_ic[temp])

        # ttstep <- 0.001
        # ttgrid <- seq(time0,tau,by=ttstep)
        # slamtA0_c_bd[it,jj,ics] <- sum(epanker(tvalue, time[temp], b)*wt[temp]/Sf_ic[temp])/sum(epanker(tvalue, ttgrid, b)*ttstep)


      }
    }
  }


  icov <- trtpos# test the treatment effect (trtpos)
  if (sum(is.na(sbeta_ic[icov,])) == 0){
    ## hypothesis testing

    alphanull <- rep(log(1-VEnull),ncs)
    alphanull.aprch2 <- alphanull[missmodel]

    alpha <- sbeta_ic[icov,]
    alpha.aprch2 <- alpha[missmodel]

    cov_ic <- covi(time,covar2,cause,sbeta_ic,wipw,strata.num,maxit,nstrt,dr,Ipsi,Spsi)

    cov_alpha <- matrix(0, ncs, ncs)
    for (ics in 1:ncs){
      for (jcs in 1:ncs){
        cov_alpha[ics,jcs] <- cov_ic[(ics-1)*ncov+icov,(jcs-1)*ncov+icov]
      }
    }
    cov_alpha.aprch2 <- cov_alpha[missmodel,missmodel]

    sigma <- sstd_ic[icov,]
    sigma.aprch2 <- sigma[missmodel]


    # approximate the distribution of test statistics
    D <- 1000
    tuta <- mvrnorm(D,rep(0,ncs),cov_alpha)
    tuta.aprch2 <- tuta[,missmodel]

    U_tuta <- tuta%*%diag(1/sigma)
    U_tuta.aprch2 <- tuta.aprch2%*%diag(1/sigma.aprch2)

    U1_tuta <- apply(U_tuta,1,FUN = min)
    U1_tuta.aprch2 <- apply(U_tuta.aprch2,1,FUN = min)
    U2_tuta <- rowSums(U_tuta^2)
    U2_tuta.aprch2 <- rowSums(U_tuta.aprch2^2)

    # H0: VE(j) <= VEnull
    # HA1 \alpha_j <= c0
    U1 <- min((alpha-alphanull)/sigma)
    pval.A1 <- mean(U1_tuta<U1)
    U1.aprch2 <- min((alpha.aprch2-alphanull.aprch2)/sigma.aprch2)
    pval.A1.aprch2 <- mean(U1_tuta.aprch2<U1.aprch2)


    # H0: VE(j) = VEnull
    # HA2 \alpha_j \ne c0
    U2 <- sum(((alpha-alphanull)/sigma)^2)
    pval.A2 <- mean(U2_tuta>U2)
    U2.aprch2 <- sum(((alpha.aprch2-alphanull.aprch2)/sigma.aprch2)^2)
    pval.A2.aprch2 <- mean(U2_tuta.aprch2>U2.aprch2)


    # H0: VE(0) <= VEnull
    # H0: VE(1) <= VEnull
    # HA1 HA2 for each cause j
    # HA1 \alpha_j <= c0
    U1j <- (alpha-alphanull)/sigma
    pval.A1j <- colMeans(U_tuta<matrix(rep(U1j,D),nrow=D,byrow=T))


    # H0: VE(0) = VEnull
    # H0: VE(1) = VEnull
    # HA2 \alpha_j \ne c0
    U2j <- ((alpha-alphanull)/sigma)^2
    pval.A2j <- colMeans(U_tuta^2>matrix(rep(U2j,D),nrow=D,byrow=T))



    # H0: VE(0) <= VE(1)
    # HB1 \alpha_1 <=...<= \alpha_J
    # w1=wJ+1=0, wj=1 for j=2,...,J-1
    w <- c(0,rep(1, ncs-1),0)
    omega <- matrix(0, ncs, ncs-1)
    for (ics in 2:ncs){
      omega[ics-1,ics-1] <- -w[ics]
      omega[ics,ics-1] <- w[ics]
    }
    diffsigma <- sqrt(diag(t(omega)%*%cov_alpha%*%omega))


    T1 <- min(alpha%*%omega/diffsigma)
    if (ncs>2){
      T_tuta <- tuta%*%omega%*%diag(1/diffsigma)
    } else {
      T_tuta <- tuta%*%omega/diffsigma
    }

    T1_tuta <- apply(T_tuta,1,FUN = min)
    pval.B1 <- mean(T1_tuta>T1)


    # H0: VE(0) = VE(1)
    # HB2 \alpha_i \ne \alpha_j, i<j
    T2 <- sum((alpha%*%omega/diffsigma)^2)
    T2_tuta <- rowSums(T_tuta^2)
    pval.B2 <- mean(T2_tuta>T2)

    ## Point estimate and 95# CI for VE
    sVE_ic <- 1-exp(alpha)
    sVEstd_ic <- sigma*exp(alpha)
    # VD=(1-VE(j)) / (1-VE(j-1))
    sVD_ic <- as.vector(exp(alpha%*%omega))
    sVDstd_ic <- as.vector(diffsigma*exp(alpha%*%omega))


    # aprch2 for HB
    ncs.aprch2 <- sum(missmodel)
    # H0: VE(0) <= VE(1)
    # HB1 \alpha_1 <=...<= \alpha_J
    # w1=wJ+1=0, wj=1 for j=2,...,J-1
    w <- c(0,rep(1, ncs.aprch2-1),0)
    omega <- matrix(0, ncs.aprch2, ncs.aprch2-1)
    for (ics in 2:ncs.aprch2){
      omega[ics-1,ics-1] <- -w[ics]
      omega[ics,ics-1] <- w[ics]
    }
    diffsigma.aprch2 <- sqrt(diag(t(omega)%*%cov_alpha.aprch2%*%omega))


    T1.aprch2 <- min(alpha.aprch2%*%omega/diffsigma.aprch2)
    if (ncs.aprch2>2){
      T_tuta.aprch2 <- tuta.aprch2%*%omega%*%diag(1/diffsigma.aprch2)
    } else {
      T_tuta.aprch2 <- tuta.aprch2%*%omega/diffsigma.aprch2
    }

    T1_tuta.aprch2 <- apply(T_tuta.aprch2,1,FUN = min)
    pval.B1.aprch2 <- mean(T1_tuta.aprch2>T1.aprch2)


    # H0: VE(0) = VE(1)
    # HB2 \alpha_i \ne \alpha_j, i<j
    T2.aprch2 <- sum((alpha.aprch2%*%omega/diffsigma.aprch2)^2)
    T2_tuta.aprch2 <- rowSums(T_tuta.aprch2^2)
    pval.B2.aprch2 <- mean(T2_tuta.aprch2>T2.aprch2)

  }

  res <- list(causes=causelevels,
              misscauses=causelevels[missmodel],
              coef.cc=sbeta_c,
              se.cc=sstd_c,
              coef=sbeta_ic,
              se=sstd_ic,
              coef.VE=sVE_ic,
              se.VE=sVEstd_ic,
              coef.VD=sVD_ic,
              se.VD=sVDstd_ic,
              U1=U1.aprch2,
              U2=U2.aprch2,
              U1j=U1j,
              U2j=U2j,
              T1=T1.aprch2,
              T2=T2.aprch2,
              pval.A1=pval.A1.aprch2,
              pval.A2=pval.A2.aprch2,
              pval.A1j=pval.A1j,
              pval.A2j=pval.A2j,
              pval.B1=pval.B1.aprch2,
              pval.B2=pval.B2.aprch2,
              tgrid=tgrid,
              lambda0=slamtA0_c,
              covariates=colnames(covar2),
              trtpos=trtpos,
              VEnull=VEnull,
              diffsigma=diffsigma,
              cov.alpha=cov_alpha)
  class(res) <- "markPH.ipw2"
  return(res)
}


##' @export
print.markPH.ipw2 <- function(x, digits=4,...){

  ncs <- ncol(x$coef)
  ncov <- nrow(x$coef)

  cat("Table 1: estimates of covaraite coefficients\n")
  for (ics in 1:ncs){
    cat("Cause:", x$causes[ics], "\n")
    table1 <- cbind(x$coef[,ics],
                    x$se[,ics],
                    pnorm(-abs(x$coef[,ics]/x$se[,ics]))*2,
                    x$coef.cc[,ics],
                    x$se.cc[,ics],
                    pnorm(-abs(x$coef.cc[,ics]/x$se.cc[,ics]))*2)
    colnames(table1) <- c("coef",
                          "se",
                          "pval",
                          "coef.cc",
                          "se.cc",
                          "pval.cc")
    rownames(table1) <- x$covariates
    print.default(table1, digits=digits)
    cat("\n")
  }

  cat("Table 2: results for VE\n")
  table2 <- cbind(x$coef.VE,x$se.VE,
                  1-(1-x$coef.VE)*exp(qnorm(0.975)*x$se[x$trtpos,]),
                  1-(1-x$coef.VE)*exp(-qnorm(0.975)*x$se[x$trtpos,]),
                  x$U1j,
                  x$pval.A1j,
                  x$U2j,
                  x$pval.A2j)
  colnames(table2) <- c("est",
                        "se",
                        "95% LL",
                        "95% UL",
                        "U1j",
                        "p-value for HAj1",
                        "U2j",
                        "p-value for HAj2")
  rownames(table2) <- paste0("VE(",x$causes,")")
  print.default(table2, digits=digits)
  cat(paste0("HAj1: VE(j)>",x$VEnull),"\n")
  cat(paste0("HAj2: VE(j) not equal to ",x$VEnull),"\n")
  cat("\n")

  cat("Table 3: results for VD\n")
  cat("VD(j,k)=(1-VE(j))/(1-VE(k))\n")
  table3 <- rbind(cbind(x$coef.VD,x$se.VD,x$coef.VD*exp(-qnorm(0.975)*x$diffsigma),x$coef.VD*exp(qnorm(0.975)*x$diffsigma)),cbind(1/x$coef.VD,x$se.VD/x$coef.VD^2,1/x$coef.VD*exp(-qnorm(0.975)*x$diffsigma),1/x$coef.VD*exp(qnorm(0.975)*x$diffsigma)))
  colnames(table3) <- c("est",
                        "se",
                        "95% LL",
                        "95% UL")
  rnames <- rep(0, (ncs-1)*2)
  for (ics in 2:ncs){
    rnames[ics-1] <- paste0("VD(",x$causes[ics],",",x$causes[ics-1],")")
    rnames[ncs-1+ics-1] <- paste0("VD(",x$causes[ics-1],",",x$causes[ics],")")
  }
  rownames(table3) <- rnames
  print.default(table3, digits=digits)
  cat("\n")

  cat("Table 4: results for hypothesis tests\n")
  cat("HA1: VE(j)>=", x$VEnull, " with strict inequality for some j in (",paste(x$misscauses,collapse = ","),")\n",sep="")
  cat("U1=",x$U1,"p-value=", x$pval.A1, "\n")
  cat("HA2: VE(j) not equal to ", x$VEnull, " for some j in (",paste(x$misscauses,collapse = ","),")\n",sep="")
  cat("U2=",x$U2,"p-value=", x$pval.A2, "\n")
  cat("HB1: VE(",head(x$misscauses,1),")>=...>=VE(",tail(x$misscauses,1),") with at least one strict inequality\n",sep="")
  cat("T1=",x$T1,"p-value=", x$pval.B1, "\n")
  cat("HB2: VE(i) not equal to VE(j) for at least one pair of i and j in (",paste(x$misscauses,collapse = ","),")\n",sep="")
  cat("T2=",x$T2,"p-value=", x$pval.B2, "\n")
  cat("\n")
}


