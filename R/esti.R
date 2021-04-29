esti <- function(time,covar,deltacs,beta,wipw,strata.num,maxit,nstrt,dr,Ipsi,Spsi){

  nsamp <- nrow(covar)
  ncov <- ncol(covar)
  npsi <- ncol(dr);

  error <- 1
  iter <- 1
  while (error>0.00001 & iter<maxit){
    U <- rep(0, ncov)
    FF <- matrix(0, ncov, ncov)
    #FU <- matrix(0, ncov, ncov)

    ebz <- exp(covar%*%beta)
    wt <- wipw*deltacs

    for (i in 1:nsamp){
      if (wt[i] > 0){

        temp <- (time >= time[i])*(strata.num==strata.num[i])
        s0 <- sum(ebz*wipw*temp)
        s1 <- rep(0,ncov)
        s2 <- matrix(0,ncov,ncov)
        for (k in 1:ncov){
          s1[k] <- sum(ebz*covar[,k]*wipw*temp)
          for (j in 1:ncov){
            s2[j,k] <- sum(ebz*covar[,j]*covar[,k]*wipw*temp)
          }
        }

        for (k in 1:ncov){
          U[k] <- U[k]+wt[i]*(covar[i,k]-s1[k]/s0)
          for (j in 1:ncov){
            FF[j,k] <- FF[j,k]+wt[i]*(s2[j,k]/s0-s1[j]*s1[k]/s0^2)
            #FU[j,k] <- FU[j,k]+(wt[i])^2*(covar[i,j]-s1[j]/s0)*(covar[i,k]-s1[k]/s0)
          }
        }
      }
    }
    #try()
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
    U2 <- matrix(0,nsamp, ncov)
    D <- array(0,c(ncov, npsi, nstrt))

    ebz <- exp(covar%*%betahat)
    wt <- wipw*deltacs

    s0 <- rep(0,nsamp)
    s1 <- matrix(0,nsamp,ncov)
    for (i in 1:nsamp){
      if (wt[i]>0){
        temp <- (time >= time[i])*(strata.num==strata.num[i])
        s0[i] <- sum(ebz*wipw*temp)
        for (k in 1:ncov){
          s1[i,k] <- sum(ebz*covar[,k]*wipw*temp)
        }
      }
    }

    for (i in 1:nsamp){
      D1 <- matrix(0,ncov,npsi)
      if (wt[i]>0){
        D1 <- -wt[i]^2*outer((covar[i,]-s1[i,]/s0[i]),dr[i,])
      }

      D2 <- matrix(0,ncov,npsi)

      for (k in 1:ncov){
        temp <- which((s0!=0)&(time<=time[i])&(strata.num==strata.num[i]))
        D2[k,] <- -wipw[i]^2*ebz[i]*sum((covar[i,k]-s1[temp,k]/s0[temp])*wt[temp]/s0[temp])*dr[i,]
      }

      D[,,strata.num[i]]=D[,,strata.num[i]]+(D1-D2)/sum(strata.num==strata.num[i])
    }

    for (i in 1:nsamp){
      U21 <- rep(0,ncov)
      if (wt[i]>0){
        U21 <- wt[i]*(covar[i,]-s1[i,]/s0[i])
      }

      U22 <- rep(0,ncov)

      for (k in 1:ncov){
        temp <- which((s0!=0)&(time<=time[i])&(strata.num==strata.num[i]))
        U22[k] <- wipw[i]*ebz[i]*sum((covar[i,k]-s1[temp,k]/s0[temp])*wt[temp]/s0[temp])
      }

      U2[i,] <- U21-U22+D[,,strata.num[i]]%*%ginv(Ipsi[,,strata.num[i]])%*%Spsi[i,]
    }
    FU <- crossprod(U2)

    # try()
    varhat <- ginv(FF)%*%FU%*%ginv(FF)
    betasig <- sqrt(diag(varhat))
  } else {
    betasig <- rep(NA,ncov)
  }

  return(list(est=betahat, se=betasig))
}









