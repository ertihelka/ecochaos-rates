# Testing GPDD time series for nonlinear stochasticity

# Load libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here, rdist, tseriesChaos, fractal)
#The fractal package has been archived now and so you may have to download
#it manually from the CRAN archive

# List all time series
here::i_am("stochasticity_tests_code.R")
input_dir <- setwd(here("data"))
input_files <- list.files(pattern = "*.csv")

# Loop and run tests on each
for(i in seq_along(input_files)){
  df <- read.csv(input_files[i])
  
  #Step 1: Select time series to analyse
  x <- as.numeric(df$PopRescale)
  x <- na.omit(x)
  
  print("**************************************")
  print(input_files[i])
  
  #Step 2: Embed time series (x)
  dump("embed_udf", file=here::here('embed_udf.R'))
  source(here::here('embed_udf.R'))
  results.embed_udf<-embed_udf(x)
  d<-results.embed_udf[[1]]
  m<-results.embed_udf[[2]]
  tw<-results.embed_udf[[3]] #embedding dimension
  Mx<-results.embed_udf[[4]] #embedded data matrix
  dist<-rdist(Mx) #distance matrix, rdist(fields)
  
  #Step 3: Nonlinear prediction skill
  source(here::here('predict_np_udf.R'))
  results.np<-predict_np_udf(Mx,frac.learn=0.5)
  nse.x<-results.np[[1]]
  print("nonlinear prediction skill:")
  print(nse.x) #nonlinear prediction skill (Nash-Sutcliffe model efficiency)
  
  #Step 4: Kaplan's delta-epsilon test for determinism
  kap.test <- determinism(x) # delta-epsilon test
  #print(kap.test) # summary of the analysis
  file_name = paste("Kaplan_test_GPDD_", tools::file_path_sans_ext(input_files[i]), 
                    ".tiff", sep="") # name plots
  tiff(file_name, units="in", width=5, height=5, res=300)
  print(plot(kap.test)) # plots E-statistics of  original series and surrogates
  dev.off()
}
