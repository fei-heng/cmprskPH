#' AIPW Estimation and Hypothesis Testing of Strain-Specific
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
#' @param missformula a formula object started with a '~' operator, and the independent
#' terms on the right of '~' are variables used for predicting the probability
#' of being not missing under a logistic regression model.
#' @param markformula a formula object started with a '~' operator, and the independent
#' terms on the right of '~' are variables used for predicting the probability
#' of being cause j under a multinomial log-linear regression model (function
#' multinom() in nnet package).
#' @param data a data.frame with the variables.
#' @param VEnull the assumed VE value in the null hypothesis. The default value
#' is 0.3.
#' @param maxit Maximum number of iterations to attempt for convergence. The
#' default is 15.
#' @param ipw Whether to conduct the IPW estimation. The default value is TRUE.
#' @param cc Whether to conduct the complete-case estimation. The default value is TRUE.
#' @param lambda0 Whether to estimate the cumulative baseline hazard. The default value is TRUE.
#'
#' @return returns an object of type 'markPH.aipw'. With the following arguments:
#' \item{causes}{the types of causes of failure}
#' \item{coef.cc}{estimates of covariate coefficients using a complete-case analysis}
#' \item{se.cc}{estimates of standard errors of estimators for covariate coefficients
#' using a complete-case analysis}
#' \item{coef.ipw}{estimates of covariate coefficients using IPW method}
#' \item{se.ipw}{estimates of standard errors of estimators for covariate coefficients
#' using IPW method}
#' \item{coef}{estimates of covariate coefficients using AIPW method}
#' \item{se}{estimates of standard errors of estimators for covariate coefficients
#' using AIPW method}
#' \item{coef.VE}{estimates of the strain-specific vaccine efficacy to reduce
#' susceptibility to strain j using AIPW method, j=1,...,J}
#' \item{se.VE}{estimates of standard errors of estimators for the strain-specific
#' vaccine efficacy to reduce susceptibility to strain j using AIPW method, j=1,...,J}
#' \item{coef.VD}{estimates of VD using AIPW method}
#' \item{se.VD}{estimates of standard errors of estimators for VD using AIPW method}
#' \item{U1}{test statistic for testing against the alternative hypothesis
#' HA1: VE(j)>=VEnull with strict inequality for some 1<=j<=J}
#' \item{U2}{test statistic for testing against the alternative hypothesis
#' HA2: VE(j) is not equal to VEnull for some 1<=j<=J}
#' \item{U1j}{test statistic for testing against the alternative hypothesis
#' HAj1: VE(j)>VEnull}
#' \item{U2j}{test statistic for testing against the alternative hypothesis
#' HAj2: VE(j) is not equal to VEnull}
#' \item{T1}{test statistic for testing against the alternative hypothesis
#' HB1: VE(1)>=...>=VE(J) with at least one strict inequality}
#' \item{T2}{test statistic for testing against the alternative hypothesis
#' HB2: VE(i) is not equal to VE(j) for at least one pair of i and j, 1<=i<j<=J}
#' \item{pval.A1}{p value for testing against the alternative hypothesis
#' HA1: VE(j)>=VEnull with strict inequality for some 1<=j<=J}
#' \item{pval.A2}{p value for testing against the alternative hypothesis
#' HA2: VE(j) is not equal to VEnull for some 1<=j<=J}
#' \item{pval.A1j}{p value for testing against the alternative hypothesis
#' HAj1: VE(j)>VEnull}
#' \item{pval.A2j}{p value for testing against the alternative hypothesis
#' HAj2: VE(j) is not equal to VEnull}
#' \item{pval.B1}{p value for testing against the alternative hypothesis
#' HB1: VE(1)>=...>=VE(J) with at least one strict inequality}
#' \item{pval.B2}{p value for testing against the alternative hypothesis
#' HB2: VE(i) is not equal to VE(j) for at least one pair of i and j, 1<=i<j<=J}
#' \item{tgrid}{specifies the time grids for estimating \eqn{\lambda_{0kj}(\cdot)}
#' , the unspecified cause-specific baseline hazard function
#' for cause j and the k-th stratum.}
#' \item{lambda0}{estimates of \eqn{\lambda_{0kj}(\cdot)} using AIPW method}
#' \item{covariates}{covariates in the cause-specific proportional hazards model.}
#' \item{trtpos}{the position of the treatment group indicator in the RHS of
#' the formula object cmprskPHformula.}
#' \item{VEnull}{the assumed VE value in the null hypothesis.}
#' \item{diffsigma}{estimates of standard deviation of \eqn{\hat\alpha_j-\hat\alpha_{j-1}}, j=2,...,J.}
#' \item{cov.alpha}{estimated covariance matrix of \eqn{(\hat\alpha_1,\dots,\hat\alpha_J)}.}
#' \item{pval.C}{p values for testing pairwise multiple comparisons of VE.}
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
#' res.aipw <- markPH.aipw(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
#'                         trtpos=1,
#'                         strata="strata",
#'                         causelevels=c(1,2),
#'                         missformula=~z1+A,
#'                         markformula=~time+z1+A,
#'                         data=sim2cs,
#'                         VEnull=0.3,
#'                         maxit=15)
#' res.aipw # print the result
#'
#' ## Example 2: Simulated competing risks data with 3 causes subject to missingness
#' data(sim3cs)
#'
#' res.aipw <- markPH.aipw(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
#'                         trtpos=1,
#'                         strata="strata",
#'                         causelevels=c(1,2,3),
#'                         missformula=~z1+A,
#'                         markformula=~time+z1+A,
#'                         data=sim3cs,
#'                         VEnull=0.3,
#'                         maxit=15)
#' res.aipw # print the result
#'
#' @import MASS nnet caret
#'
#' @export




markPH.aipw <- function(cmprskPHformula,
                        trtpos=1,
                        strata,
                        causelevels=NULL,
                        missformula,
                        markformula,
                        data=parent.frame(),
                        VEnull=0.3,
                        maxit=15,
                        ipw=T,
                        cc=T,
                        lambda0=T){
  # next version:
  # ...

  data[data=="NaN"] <- NA


  a <- model.frame(cmprskPHformula, data=data, na.action = na.pass)
  time <- model.response(a)[,1]
  delta <- model.response(a)[,2] # 0: censored, 1: observed failure
  if (is.null(causelevels)){
    cause.fa <- factor(model.response(a)[,3])
  } else {
    # define causes following the order of causelevels
    cause.fa <- factor(model.response(a)[,3], labels=causelevels)
  }
  causelevels <- levels(cause.fa)
  cause <- as.numeric(cause.fa) # cause = 1,2,...,J
  cause[is.na(cause.fa)] <- NA
  data$cause <- cause

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
  markformula <- formula(paste(c("cause",markformula),collapse = " "))
  ## predict P(R=1) and P(cause=k) in jj-th stratum
  # predict P(R=1) (in R, use subset=)
  # logistic regression
  wipw <- R
  a <- model.frame(missformula, data=data, na.action = na.pass)
  covar.miss <- as.matrix(model.matrix(a, data=a, na.action = na.pass)[,-1])
  colnames(covar.miss) <- colnames(a)[-1]

  npsi <- ncol(covar.miss)+1
  dr <- matrix(0,nsamp,npsi)
  Ipsi <- array(0, c(npsi,npsi,nstrt))
  Spsi <- matrix(0,nsamp,npsi)
  for (jj in 1:nstrt){

    temp <- which((delta==1)&(strata.num==jj))
    miss.res <- glm(missformula, data=data, family='binomial',subset=temp)
    newdata <- as.data.frame(covar.miss[temp,])
    colnames(newdata) <- colnames(a)[-1]
    rhat <- predict(miss.res, newdata, type="response")

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

  # predict cause/mark V for delta==1
  # multinomial logistic regression
  # library(nnet)
  a <- model.frame(markformula, data=data, na.action = na.pass)
  covar.mark <- as.matrix(model.matrix(a, data=a, na.action = na.pass)[,-1])

  rhohat <- matrix(0,nsamp,ncs)
  colnames(rhohat) <- causelevels
  if (ncs==2){
    for (jj in 1:nstrt){


      # ## method 1
      # # multinom() from package nnet
      # # Fits multinomial log-linear models via neural networks
      # # other options: mlogit() from package mlogit
      # temp <- (delta==1)&(strata.num==jj)&R
      # mark.res <- multinom(markformula, data=data, subset=temp, trace=F)
      # # test for ncs>2!!!
      # # newdata <- as.data.frame(covar.mark[(delta==1)&(strata.num==jj),])
      # # colnames(newdata) <- colnames(a)[-1]
      # newdata <- data[(delta==1)&(strata.num==jj), ]
      # rhohat[(delta==1)&(strata.num==jj),2] <- predict(mark.res, newdata, type="probs")

      # method 2 knn

      proba <- predict(knn3(covar.mark[(delta==1)&(strata.num==jj)&R,],
                            factor(data$cause[(delta==1)&(strata.num==jj)&R]),
                            k = 5),
                       covar.mark[(delta==1)&(strata.num==jj),])
      proba.class <- colnames(rhohat)[as.numeric(colnames(proba))]
      rhohat[(delta==1)&(strata.num==jj),proba.class] <- proba
    }
  } else {
    for (jj in 1:nstrt){
      # method 1
      # temp <- (delta==1)&(strata.num==jj)&R
      # # multinom() from package nnet
      # # Fits multinomial log-linear models via neural networks
      # # other options: mlogit() from package mlogit
      # mark.res <- multinom(markformula, data=data, subset=temp, trace=F)
      # # test for ncs>2!!!
      # newdata <- as.data.frame(covar.mark[(delta==1)&(strata.num==jj),])
      # colnames(newdata) <- colnames(a)[-1]
      # rhohat[(delta==1)&(strata.num==jj), mark.res$lab] <- predict(mark.res, newdata, type="probs")

      # method 2 knn

      proba <- predict(knn3(covar.mark[(delta==1)&(strata.num==jj)&R,],
                            factor(data$cause[(delta==1)&(strata.num==jj)&R]),
                            k = 5),
                       covar.mark[(delta==1)&(strata.num==jj),])
      proba.class <- colnames(rhohat)[as.numeric(colnames(proba))]
      rhohat[(delta==1)&(strata.num==jj),proba.class] <- proba
    }
  }


  ## estimation for each cause
  # initialization
  sbeta_c <- matrix(NA, ncov, ncs)
  sstd_c <- matrix(NA, ncov, ncs)

  sbeta_ic <- matrix(NA, ncov, ncs)
  sstd_ic <- matrix(NA, ncov, ncs)

  sbeta_acc <- matrix(NA, ncov, ncs)
  sstd_acc <- matrix(NA, ncov, ncs)


  sVE_acc <- rep(NA,ncs)
  sVEstd_acc <- rep(NA,ncs)
  sVD_acc <- rep(NA,ncs-1)
  sVDstd_acc <- rep(NA,ncs-1)


  for (ics in 1:ncs){
    # deltacs: indicator of non-censoring for each cause
    deltacs <- cause == ics# be carefule: NaN ->F
    deltacs[is.na(cause)] <- F

    # complete estimation
    if (cc==T){
      beta0 <- rep(0,ncov)
      res.cc <- estf(time,covar2,deltacs,beta0,strata.num,maxit,subset=R)
      # complete result
      sbeta_c[,ics] <- res.cc$est
      sstd_c[,ics] <- res.cc$se
    }


    # ipw estimation
    if (ipw==T){
      beta0 <- rep(0,ncov)
      res.ipw <- esti(time,covar2,deltacs,beta0,wipw,strata.num,maxit,nstrt,dr,Ipsi,Spsi)
      # ipw-c result
      sbeta_ic[,ics] <- res.ipw$est
      sstd_ic[,ics] <- res.ipw$se
    }


    # aipw estimation
    beta0 <- rep(0,ncov)
    res.aipw <- esta(time,covar2,deltacs,beta0,wipw,rhohat[,ics],delta,strata.num,maxit)

    # aipw_cc result
    sbeta_acc[,ics] <- res.aipw$est
    sstd_acc[,ics] <- res.aipw$se
  }

  ## calculate the cumulative baseline function
  # aipw
  # sLamtA0_c <- array(0, c(ngrid,nstrt,ncs))
  slamtA0_c <- array(0, c(ngrid,nstrt,ncs))
  # slamtA0_c_bd <- array(0, c(ngrid,nstrt,ncs)) boundary correction
  if (lambda0==T){
    for (ics in 1:ncs){
      deltacs <- cause == ics# be carefule: NaN ->F
      deltacs[is.na(cause)] <- F


      ebz <- exp(covar2%*%sbeta_acc[,ics])

      Sf_acc <- rep(0, nsamp)
      for (i in 1:nsamp){
        if (delta[i]){
          Sf_acc[i] <- sum(ebz*(time >= time[i])*(strata.num==strata.num[i]))
        }
      }

      for (jj in 1:nstrt){
        for (it in 1:ngrid){
          tvalue <- tgrid[it]
          wt <- (wipw*deltacs+(1-wipw)*rhohat[,ics])*delta

          # temp=which((wt!=0)&(Sf_acc!=0)&(time<tvalue)&(strata.num==jj))
          # sLamtA0_c[it,jj,ics]  <- sum(wt[temp]/Sf_acc[temp])

          temp=which((wt!=0)&(Sf_acc!=0)&(strata.num==jj))
          slamtA0_c[it,jj,ics] <- sum(epanker(tvalue, time[temp], b)*wt[temp]/Sf_acc[temp])

          # ttstep <- 0.001
          # ttgrid <- seq(time0,tau,by=ttstep)
          # slamtA0_c_bd[it,jj,ics] <- sum(epanker(tvalue, time[temp], b)*wt[temp]/Sf_acc[temp])/sum(epanker(tvalue, ttgrid, b)*ttstep)


        }
      }
    }
  }



  ## hypothesis tests for causes with non-NA estimates
  icov <- trtpos# test the treatment effect (trtpos)
  U1 <- U2 <- T1 <- T2 <- pval.A1 <- pval.A2 <- pval.B1 <- pval.B2 <- NA
  U1j <- U2j <- pval.A1j <- pval.A2j <- rep(NA, ncs)
  diffsigma <- rep(NA, ncs)


  # causes <- levels(cause.fa)
  causes <- 1:ncs
  nonna <- !is.na(sbeta_acc[icov,])
  causes.na <- causes[!nonna]

  id.test <- !(cause%in%causes.na)
  time.test <- time[id.test]
  covar2.test <- covar2[id.test,,drop=F]
  cause.test <- cause[id.test]
  sbeta_acc.test <- sbeta_acc[,nonna,drop=F]
  wipw.test <- wipw[id.test]
  rhohat.test <- rhohat[id.test,nonna]
  delta.test <- delta[id.test]
  strata.num.test <- strata.num[id.test]

  ncs.test <- sum(nonna)
    alphanull <- rep(log(1-VEnull),ncs.test)

    cov_acc <- cova(time.test,covar2.test,cause.test,sbeta_acc.test,wipw.test,rhohat.test,delta.test,strata.num.test,maxit,causes[nonna])

    alpha <- sbeta_acc[icov,nonna]
    cov_alpha <- matrix(NA, ncs.test, ncs.test)
    for (ics in 1:ncs.test){
      for (jcs in 1:ncs.test){
        cov_alpha[ics,jcs] <- cov_acc[(ics-1)*ncov+icov,(jcs-1)*ncov+icov]
      }
    }

    sigma <- sstd_acc[icov,nonna]

    # approximate the distribution of test statistics
    D <- 1000
    tuta <- mvrnorm(D,rep(0,ncs.test),cov_alpha)
    U_tuta <- tuta%*%diag(1/sigma)
    U1_tuta <- apply(U_tuta,1,FUN = min)
    U2_tuta <- rowSums(U_tuta^2)

    # H0: VE(j) <= VEnull
    # HA1 \alpha_j <= c0
    U1 <- min((alpha-alphanull)/sigma)
    pval.A1 <- mean(U1_tuta<U1)


    # H0: VE(j) = VEnull
    # HA2 \alpha_j \ne c0
    U2 <- sum(((alpha-alphanull)/sigma)^2)
    pval.A2 <- mean(U2_tuta>U2)


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
    w <- c(0,rep(1, ncs.test-1),0)
    omega <- matrix(0, ncs.test, ncs.test-1)
    for (ics in 2:ncs.test){
      omega[ics-1,ics-1] <- -w[ics]
      omega[ics,ics-1] <- w[ics]
    }
    diffsigma <- sqrt(diag(t(omega)%*%cov_alpha%*%omega))

    T1 <- min(alpha%*%omega/diffsigma)
    if (ncs.test>2){
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
    sVE_acc <- 1-exp(alpha)
    sVEstd_acc <- sigma*exp(alpha)
    # VD=(1-VE(j)) / (1-VE(j-1))
    sVD_acc <- as.vector(exp(alpha%*%omega))
    sVDstd_acc <- as.vector(diffsigma*exp(alpha%*%omega))

    ## added by Fei - 5/22/2023
    # H0: VE(i) = VE(j)
    # HC: VE(i) \ne VE(j)
    pval.C <- matrix(0, nrow=ncs.test, ncol=ncs.test)
    for (ics in 1:(ncs.test-1)){
      for(jcs in (ics+1):ncs.test){
        T2ij <- (alpha[ics]-alpha[jcs])^2/(cov_alpha[ics,ics]-2*cov_alpha[ics,jcs]+cov_alpha[jcs,jcs])
        T2ij_tuta <- (tuta[,ics]-tuta[,jcs])^2/(cov_alpha[ics,ics]-2*cov_alpha[ics,jcs]+cov_alpha[jcs,jcs])
        pval.C[ics, jcs] <- mean(T2ij_tuta>T2ij)
      }
    }



  res <- list(causes=causes[nonna],
              causes.na=causes.na,
              coef.cc=sbeta_c[,nonna,drop=F],
              se.cc=sstd_c[,nonna,drop=F],
              coef.ipw=sbeta_ic[,nonna,drop=F],
              se.ipw=sstd_ic[,nonna,drop=F],
              coef=sbeta_acc[,nonna,drop=F],
              se=sstd_acc[,nonna,drop=F],
              coef.VE=sVE_acc,
              se.VE=sVEstd_acc,
              coef.VD=sVD_acc,
              se.VD=sVDstd_acc,
              U1=U1,
              U2=U2,
              U1j=U1j,
              U2j=U2j,
              T1=T1,
              T2=T2,
              pval.A1=pval.A1,
              pval.A2=pval.A2,
              pval.A1j=pval.A1j,
              pval.A2j=pval.A2j,
              pval.B1=pval.B1,
              pval.B2=pval.B2,
              tgrid=tgrid,
              lambda0=slamtA0_c[,,nonna,drop=F],
              covariates=colnames(covar2),
              trtpos=trtpos,
              VEnull=VEnull,
              diffsigma=diffsigma,
              cov.alpha=cov_alpha,
              pval.C=pval.C)
  class(res) <- "markPH.aipw"
  return(res)
}


##' @export
print.markPH.aipw <- function(x, digits=4,...){

  ncs <- ncol(x$coef)
  ncov <- nrow(x$coef)

  cat("\n*********************************************************************\n")
  cat("Causes {", x$causes.na, "} with NA estimates are ignored in the results\n")
  cat("*********************************************************************\n\n")

  cat("Table 1: estimates of covaraite coefficients\n")
  for (ics in 1:ncs){
    cat("Cause:", x$causes[ics], "\n")
    table1 <- cbind(x$coef[,ics],
                    x$se[,ics],
                    pnorm(-abs(x$coef[,ics]/x$se[,ics]))*2,
                    x$coef.ipw[,ics],
                    x$se.ipw[,ics],
                    pnorm(-abs(x$coef.ipw[,ics]/x$se.ipw[,ics]))*2,
                    x$coef.cc[,ics],
                    x$se.cc[,ics],
                    pnorm(-abs(x$coef.cc[,ics]/x$se.cc[,ics]))*2)
    colnames(table1) <- c("coef",
                          "se",
                          "pval",
                          "coef.ipw",
                          "se.ipw",
                          "pval.ipw",
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
  cat("HA1: VE(j)>=", x$VEnull, "with strict inequality for some j\n")
  cat("U1=",x$U1,"p-value=", x$pval.A1, "\n")
  cat("HA2: VE(j) not equal to", x$VEnull, "for some j\n")
  cat("U2=",x$U2,"p-value=", x$pval.A2, "\n")
  cat("HB1: VE(1)>=...>=VE(J) with at least one strict inequality\n")
  cat("T1=",x$T1,"p-value=", x$pval.B1, "\n")
  cat("HB2: VE(i) not equal to VE(j) for at least one pair of i and j, 1<=i<j<=J\n")
  cat("T2=",x$T2,"p-value=", x$pval.B2, "\n")
  cat("\n")
}


