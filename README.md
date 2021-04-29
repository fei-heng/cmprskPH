# Cause-specific proportional hazards model for compering risks analysis (cmprskPH)

Authors: Fei Heng (f.heng@unf.edu), Yanqing Sun (yasun@uncc.edu), and Peter B. Gilbert 

R package for the paper `Estimation and Hypothesis Testing of Strain-Specific Vaccine Efficacy Adjusted for Covariate Effects with Missing Causes'

## Installation
Requires [Rtools](http://cran.r-project.org/bin/windows/Rtools/) 
on windows and [development tools](http://cran.r-project.org/bin/macosx/tools/) (+Xcode) for Mac OS X:
```{r eval=F}
devtools::install_github("fei-heng/cmprskPH")
```

## Example Codes
### Example 1: for 2 causes
```{r eval=F}
library(cmprskPH)
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
```

### Example 2: for 3 causes
```{r eval=F}
library(cmprskPH)
data(sim3cs)
res.aipw <- markPH.aipw(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
				    trtpos=1,
				    strata="strata",
				    missformula=~z1+A,
				    markformula=~time+z1+A,
				    data=sim3cs,
				    VEnull=0.3,
				    maxit=15)
res.aipw # print the result

res.ipw <- markPH.ipw(cmprskPHformula=cbind(time,delta,cause)~z1+z2,
				   trtpos=1,
				   strata="strata",
				   missformula=~z1+A,
				   data=sim3cs,
				   VEnull=0.3,
				   maxit=15)
res.ipw # print the result
```
