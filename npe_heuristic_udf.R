# A user-defined heuristic function to calculate nonlinear prediction skill
#
# Time series embedding can be computationally expensive and, particularly in R,
# may crash the session or lead to an out-of-memory error, even for short time series: 
# https://stackoverflow.com/questions/72533418/error-cannot-allocate-memory-block-of-size-67108864-tb-in-the-r-function-false
# 
# This heuristic algorithm estimates the nonlinear prediction skill of a series
# for four embedding dimensions (d = 2, 5, 10, 15) and takes a mean. This sacrifices
# accuracy for scalability to potentially large datasets, but in our experience 
# leads to satisfactory estimates, as the prediction skill is generally comparable
# across similar embedding dimensions.
#
# Erik Tihelka, 2022

npe_heuristic <- function(x) {
  source(here::here('predict_np_udf.R'))
  x <- unlist(x)
  x <- na.omit(x)
  
  #Embedding dimension: 2
  Mx <- embed(x, dimension = 2)
  source(here::here('predict_np_udf.R'))
  results.np<-predict_np_udf(Mx,frac.learn=0.5)
  npe02<-results.np[[1]]
  
  #Embedding dimension: 5
  Mx <- embed(x, dimension = 5)
  results.np<-predict_np_udf(Mx,frac.learn=0.5)
  npe05<-results.np[[1]]
  
  #Embedding dimension: 10
  Mx <- embed(x, dimension = 10)
  results.np<-predict_np_udf(Mx,frac.learn=0.5)
  npe10<-results.np[[1]]
  
  #Embedding dimension: 15
  Mx <- embed(x, dimension = 15)
  results.np<-predict_np_udf(Mx,frac.learn=0.5)
  npe15<-results.np[[1]]
  
  # Mean across the five embedding dimensions
  npe_avr <- mean(na.omit(npe02, npe05, npe10, npe15))
}