#' Estimation and Hypothesis Testing of Strain-Specific
#' Vaccine Efficacy Adjusted for Covariate Effects
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
#' @param data a data.frame with the variables.
#' @param VEnull the assumed VE value in the null hypothesis. The default value
#' is 0.3.
#' @param maxit Maximum number of iterations to attempt for convergence. The
#' default is 15.
#'
#' @return returns an object of type 'markPH.aipw'. With the following arguments:
#' \item{causes}{the types of causes of failure}
#' \item{coef}{estimates of covariate coefficients}
#' \item{se}{estimates of standard errors of estimators for covariate coefficients}
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
#'
#' @author Fei Heng
#'
#' @references (2021+)
#'
#' @examples
#'
#' ## Example 1: Simulated competing risks data with 2 causes
#' data(sim2cs)
#'
#' # complete case analysis
#' sim2cs <- sim2cs[sim2cs$R==1,]
#' sim2cs <- subset(sim2cs, select=-c(A, R))
#'
#' res <- markPH(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
#'                         trtpos=1,
#'                         strata="strata",
#'                         data=sim2cs,
#'                         VEnull=0.3,
#'                         maxit=15)
#' res # print the result
#'
#' ## Example 2: Simulated competing risks data with 3 causes
#' data(sim3cs)
#'
#' # complete case analysis
#' sim3cs <- sim3cs[sim3cs$R==1,]
#' sim3cs <- subset(sim3cs, select=-c(A, R))
#'
#' res <- markPH(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
#'                         trtpos=1,
#'                         strata="strata",
#'                         data=sim3cs,
#'                         VEnull=0.3,
#'                         maxit=15)
#' res # print the result
#'
#' @import MASS nnet
#'
#' @export




markPH <- function(cmprskPHformula,
                        trtpos=1,
                        strata,
                        data=parent.frame(),
                        VEnull=0.3,
                        maxit=15){
  # next version:
  # ...

  data[data=="NaN"] <- NA


  a <- model.frame(cmprskPHformula, data=data, na.action = na.pass)
  time <- model.response(a)[,1]
  delta <- model.response(a)[,2] # 0: censored, 1: observed failure
  cause.fa <- factor(model.response(a)[,3])
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


  ## estimation for each cause
  # initialization
  sbeta_c <- matrix(0, ncov, ncs)
  sstd_c <- matrix(0, ncov, ncs)


  for (ics in 1:ncs){
    # deltacs: indicator of non-censoring for each cause
    deltacs <- cause == ics# be carefule: NaN ->F
    deltacs[is.na(cause)] <- F

    # complete estimation
    beta0 <- rep(0,ncov)
    # library(MASS) # ginv()
    res.cc <- estf(time,covar2,deltacs,beta0,strata.num,maxit,subset=R)


    # complete result
    sbeta_c[,ics] <- res.cc$est
    sstd_c[,ics] <- res.cc$se
  }

  ## calculate the cumulative baseline function
  slamtA0_c <- array(0, c(ngrid,nstrt,ncs))
  for (ics in 1:ncs){
    deltacs <- cause == ics# be carefule: NaN ->F
    deltacs[is.na(cause)] <- F


    ebz <- exp(covar2%*%sbeta_c[,ics])

    Sf_c <- rep(0, nsamp)
    for (i in 1:nsamp){
      if (delta[i]){
        Sf_c[i] <- sum(ebz*(time >= time[i])*(strata.num==strata.num[i]))
      }
    }

    for (jj in 1:nstrt){
      for (it in 1:ngrid){
        tvalue <- tgrid[it]
        wt <- deltacs*delta

        temp=which((wt!=0)&(Sf_c!=0)&(strata.num==jj))
        slamtA0_c[it,jj,ics] <- sum(epanker(tvalue, time[temp], b)*wt[temp]/Sf_c[temp])

      }
    }
  }


  icov <- trtpos# test the treatment effect (trtpos)
  if (sum(is.na(sbeta_c[icov,])) == 0){
    ## hypothesis testing

    alphanull <- rep(log(1-VEnull),ncs)

    cov_c <- covf(time,covar2,cause,sbeta_c,
                  # wipw,rhohat, # delete
                  delta,strata.num,maxit,1:ncs)

    alpha <- sbeta_c[icov,]
    cov_alpha <- matrix(0, ncs, ncs)
    for (ics in 1:ncs){
      for (jcs in 1:ncs){
        cov_alpha[ics,jcs] <- cov_c[(ics-1)*ncov+icov,(jcs-1)*ncov+icov]
      }
    }
    # round(cov_alpha, 4)
    # Since alpha_hats are independent with each other, we can simply set cov_alpha as:
    # diag(sstd_c[icov,]^2)

    sigma <- sstd_c[icov,]

    # approximate the distribution of test statistics
    D <- 1000
    tuta <- mvrnorm(D,rep(0,ncs),cov_alpha)
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
    sVE_c <- 1-exp(alpha)
    sVEstd_c <- sigma*exp(alpha)
    # VD=(1 – VE(j)) / (1 – VE(j-1))
    sVD_c <- as.vector(exp(alpha%*%omega))
    sVDstd_c <- as.vector(diffsigma*exp(alpha%*%omega))
  }

  res <- list(causes=levels(cause.fa),
              coef=sbeta_c,
              se=sstd_c,
              coef.VE=sVE_c,
              se.VE=sVEstd_c,
              coef.VD=sVD_c,
              se.VD=sVDstd_c,
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
              lambda0=slamtA0_c,
              covariates=colnames(covar2),
              trtpos=trtpos,
              VEnull=VEnull,
              diffsigma=diffsigma,
              cov.alpha=cov_alpha)
  class(res) <- "markPH"
  return(res)
}


##' @export
print.markPH <- function(x, digits=4,...){

  ncs <- ncol(x$coef)
  ncov <- nrow(x$coef)

  cat("Table 1: estimates of covaraite coefficients\n")
  for (ics in 1:ncs){
    cat("Cause:", x$causes[ics], "\n")
    table1 <- cbind(x$coef[,ics],
                    x$se[,ics],
                    pnorm(-abs(x$coef[,ics]/x$se[,ics]))*2)
    colnames(table1) <- c("coef",
                          "se",
                          "pval")
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


