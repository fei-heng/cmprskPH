epanker <- function(tk, tvalue, hband){
  u <- (tk-tvalue)/hband
  return(3/4/hband*(1-u^2)*(abs(u)<1))
}



