# Cause-specific proportional hazards model for competing risks analysis (cmprskPH)

Authors: Fei Heng (f.heng@unf.edu), Yanqing Sun (yasun@uncc.edu), and Peter B. Gilbert 

R package for the paper `Estimation and Hypothesis Testing of Strain-Specific Vaccine Efficacy Adjusted for Covariate Effects with Missing Causes'

## Installation
Requires [Rtools](http://cran.r-project.org/bin/windows/Rtools/) 
on windows and [development tools](http://cran.r-project.org/bin/macosx/tools/) (+Xcode) for Mac OS X:
```{r eval=F}
devtools::install_github("fei-heng/cmprskPH")
```

## Code Examples
### Example 1: for 2 causes
```{r eval=F}
library(cmprskPH)

## Approach 1
data(sim2cs)
res.aipw <- markPH.aipw(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
            trtpos=1,
            strata="strata",
            missformula=~z1+A,
            markformula=~time+z1+A,
            data=sim2cs,
            VEnull=0.3,
            maxit=15)
res.aipw # print the result

res.ipw <- markPH.ipw(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
				   trtpos=1,
				   strata="strata",
				   missformula=~z1+A,
				   data=sim2cs,
				   VEnull=0.3,
				   maxit=15)
res.ipw # print the result

## Approach 2
data(sim2cs)
threshold <- 0.5
sim2cs$cause[(sim2cs$delta==1)&(sim2cs$A<threshold)] <- 3

res.aipw2 <- markPH.aipw2(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
            trtpos=1,
            strata="strata",
            causelevels=c(1,2,3),
            missmodel=c(T,T,F),
            missformula=~z1+A,
            markformula=~time+z1+A,
            data=sim2cs,
            VEnull=0.3,
            maxit=15)
res.aipw2

res.ipw2 <- markPH.ipw2(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
                   trtpos=1,
                   strata="strata",
                   causelevels=c(1,2,3),
                   missmodel=c(T,T,F),
                   missformula=~z1+A,
                   data=sim2cs,
                   VEnull=0.3,
                   maxit=15)
res.ipw2
```

### Example 2: for 3 causes
```{r eval=F}
library(cmprskPH)

## Approach 1
data(sim3cs)
res.aipw <- markPH.aipw(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
				    trtpos=1,
				    strata="strata",
				    missformula=~z1+A,
				    markformula=~time+z1+A,
				    data=sim3cs,
				    VEnull=0.3,
				    maxit=15)
res.aipw

res.ipw <- markPH.ipw(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
				   trtpos=1,
				   strata="strata",
				   missformula=~z1+A,
				   data=sim3cs,
				   VEnull=0.3,
				   maxit=15)
res.ipw

## Approach 2
data(sim3cs)
threshold <- 0.5
simdata$cause[(simdata$delta==1)&(simdata$A<threshold)] <- 4

res.aipw2 <- markPH.aipw2(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
                   trtpos=1,
                   strata="strata",
                   causelevels=c(1,2,3,4),
                   missmodel=c(T,T,T,F),
                   missformula=~z1+A,
                   markformula=~time+z1+A,
                   data=simdata,
                   VEnull=0.3,
                   maxit=15)
res.aipw2

res.ipw2 <- markPH.ipw2(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
                   trtpos=1,
                   strata="strata",
                   causelevels=c(1,2,3,4),
                   missmodel=c(T,T,T,F),
                   missformula=~z1+A,
                   data=simdata,
                   VEnull=0.3,
                   maxit=15)
res.ipw2
```
