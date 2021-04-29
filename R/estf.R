estf <- function(time,covar,deltacs,beta,strata.num,maxit,subset=NULL){
  # Be careful! This simple version only works for time-independent covariate Z

  if (!is.null(subset)){
    time <- time[subset]
    covar <- covar[subset,]
    deltacs <- deltacs[subset]
    strata.num <- strata.num[subset]
  }


  nsamp <- nrow(covar)
  ncov <- ncol(covar)

  error <- 1
  iter <- 1
  while (error>0.00001 & iter<maxit){
    U <- rep(0, ncov)
    FF <- matrix(0, ncov, ncov)

    ebz <- exp(covar%*%beta)

    wt <- deltacs

    for (i in 1:nsamp){
      if (wt[i]){

        temp <- (time >= time[i])*(strata.num==strata.num[i])
        s0 <- sum(ebz*temp)
        s1 <- rep(0,ncov)
        s2 <- matrix(0,ncov,ncov)
        for (k in 1:ncov){
          s1[k] <- sum(ebz*covar[,k]*temp)
          for (j in 1:ncov){
            s2[j,k] <- sum(ebz*covar[,j]*covar[,k]*temp)
          }
        }


        for (k in 1:ncov){
          U[k] <- U[k]+wt[i]*(covar[i,k]-s1[k]/s0)
          for (j in 1:ncov){
            FF[j,k] <- FF[j,k]+wt[i]*(s2[j,k]/s0-s1[j]*s1[k]/s0^2)
          }
        }
      }
    }

    # try()
    change <- ginv(FF)%*%U
    beta <- beta + change
    error <- sum(abs(change))
    betahat <- beta

    iter <- iter + 1
  }


  if (iter >= maxit){
    betahat <- rep(NA,ncov)
  }

  if (sum(is.na(betahat)) == 0){
    # try()
    varhat <- ginv(FF)
    betasig <- sqrt(diag(varhat))
  } else {
    betasig <- rep(NA,ncov)
  }

  return(list(est=betahat, se=betasig))
}









