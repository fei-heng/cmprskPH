cova <- function(time,covar,cause,betaall,wipw,rhohatall,delta,strata.num,maxit){

  nsamp <- nrow(covar)
  ncov <- ncol(covar)
  ncs <- length(unique(cause[!is.na(cause)]))

  U2 <- matrix(0, nsamp, ncov*ncs)
  FF <- matrix(0, ncov*ncs, ncov*ncs)


  wt <- delta

  for (ics in 1:ncs){
    pos <- ncov*(ics-1)

    deltacs <- cause == ics
    deltacs[is.na(cause)] <- F
    rhohat <- rhohatall[,ics]
    beta <- betaall[,ics]


    #D <- array(0,c(ncov, npsi, nstrt))

    ebz <- exp(covar%*%beta)
    #wt <- delta

    s0 <- rep(0,nsamp)
    s1 <- matrix(0,nsamp,ncov)
    for (i in 1:nsamp){
      if (wt[i]>0){
        temp <- (time >= time[i])*(strata.num==strata.num[i])
        s0[i] <- sum(ebz*temp)
        s2 <- matrix(0,ncov,ncov)
        for (k in 1:ncov){
          s1[i,k] <- sum(ebz*covar[,k]*temp)
          for (j in 1:ncov){
            s2[j,k] <- sum(ebz*covar[,j]*covar[,k]*temp)
          }
        }

        for (k in 1:ncov){
          for (j in 1:ncov){
            FF[pos+j,pos+k] <- FF[pos+j,pos+k]+wt[i]*(s2[j,k]/s0[i]-s1[i,j]*s1[i,k]/s0[i]^2)*(wipw[i]*deltacs[i]+(1-wipw[i])*rhohat[i])
          }
        }
      }
    }

    # for (i in 1:nsamp){
    #   D1 <- matrix(0,ncov,npsi)
    #   if (wt[i]>0){
    #     D1 <- -wt[i]^2*outer((covar[i,]-s1[i,]/s0[i]),dr[i,])
    #   }
    #
    #   D2 <- matrix(0,ncov,npsi)
    #
    #   for (k in 1:ncov){
    #     temp <- which((s0!=0)&(time<=time[i])&(strata.num==strata.num[i]))
    #     D2[k,] <- -wipw[i]^2*ebz[i]*sum((covar[i,k]-s1[temp,k]/s0[temp])*wt[temp]/s0[temp])*dr[i,]
    #   }
    #
    #   D[,,strata.num[i]]=D[,,strata.num[i]]+(D1-D2)/sum(strata.num==strata.num[i])
    # }

    for (i in 1:nsamp){
      U21 <- rep(0,ncov)
      if (wt[i]>0){
        U21 <- wt[i]*(covar[i,]-s1[i,]/s0[i])*(wipw[i]*deltacs[i]+(1-wipw[i])*rhohat[i])
      }

      U22 <- rep(0,ncov)

      for (k in 1:ncov){
        temp <- which((s0!=0)&(time<=time[i])&(strata.num==strata.num[i]))
        U22[k] <- ebz[i]*sum((covar[i,k]-s1[temp,k]/s0[temp])*(wipw[temp]*deltacs[temp]+(1-wipw[temp])*rhohat[temp])*wt[temp]/s0[temp])
      }

      U2[i,(pos+1):(pos+ncov)] <- U21-U22#+D[,,strata.num[i]]%*%ginv(Ipsi[,,strata.num[i]])%*%Spsi[i,]
    }
  }

    FU <- crossprod(U2)

  # try
  return(ginv(FF)%*%FU%*%ginv(FF))
}

