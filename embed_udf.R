# A user-defined time series embedding function
# After Huffaker, Bittelli & Rosa (2017) Nonlinear Time Series Analysis with R, 
# Oxford University Press, pp. 154-174.

embed_udf <-
function(x){
  #Step 1: Calculate embedding delay
  #Step 1a. Compute average mutual information (AMI) function
  library(tseriesChaos)
  mutual.out<-mutual(x) #mutual(tseriesChaos)
  #Step 1b: Use embedding approach to calculate delay at which AMI hits
  #its first minimum
  mutual.em<-embedd(mutual.out,2,1) #embedd(tseriesChaos)
  mutual.adj.length<-length(mutual.out)-1 #lose 1 observation to delay
  mutual.hold<-matrix(0,mutual.adj.length,1)
  for (i in 1:mutual.adj.length){ #loop to compute successfive value differences
    mutual.test<-if(mutual.em[i,1]>mutual.em[i,2])TRUE else break
    mutual.hold[i,1]<-mutual.test #"TRUE' = 1 in R, so min delay occurs at sum(TRUE)
  } #end iloop
  mutual.hold.sum<-sum(mutual.hold)
  #Step 1c: Estimate embeddeding delay (d). If the mutual information function is
  #decreasing across all 20 delays, set delay at max = 20
  d<-if(mutual.hold.sum<20)mutual.hold.sum else 20 #Embedding delay
  #Step 2: Compute Theiler window parameter (tw) required for false nearest
  # neighbours test
  #Step 2a: Autocorrelation function
  lag.max=100
  acf.run<-acf(x,lag.max)
  acf.out<-acf.run$acf #array of acf values 
  #Step 2b: Use embedding approach to calculate delay at which AMI hits 
  #its first minimum.
  acf.em<-embedd(acf.out,2,1) #embedd(tseriesChaos) 
  acf.adj.length<-length(acf.out)-1 #lose 1 observation to lag 
  acf.hold<-matrix(0,acf.adj.length,1)  
  for (i in 1:acf.adj.length){ #loop to compute successfive value differences
  acf.test<-if(acf.em[i,1]>acf.em[i,2])TRUE else break
  acf.hold[i,1]<-acf.test #"TRUE' = 1 in R, so min delay occurs at sum(TRUE)
} #end iloop
#Step 2c: Estimate Theiler window (tw)
tw<-sum(acf.hold)
#Step 3: Embedding dimension (m)
#Step 3a: False nearest neighbours function
m.max<-6 #maximum number of embedding dimensions to consider
fn.out<-false.nearest(x,m.max,d,tw) #false.nearest(tseriesChaos)
fn.out[is.na(fn.out)] <- 0 #set NA in fn.out to zero
#plot(fn.out)
#Step 3b: Find delay at which false nearest neighbours decrease below set tolerance
#Output vector of fnn percentages from fn.out
fnp<-c(fn.out[1],fn.out[3],fn.out[5],fn.out[7],fn.out[9],fn.out[11])
fnp.tol<-fnp>0.15 #If fnp greater than tolerance of 15%, T entered into fnp.tol
fnp.tol.sum<-sum(fnp.tol) #sum up number of T's
m<-if(fnp.tol.sum<m.max)fnp.tol.sum+1 else m.max #Embedding dimension
#Step 4: Embed time series (Mx)
#If m=1, embedd routine crashes due to 'subscript out of bounds' error--need to
#guarantee an embedding dimension of at least two:
if(m<=1){m<-2} else {m}
Mx<-embedd(x,m,d) #embedd(tseriesChaos)
#Results
results.embed_udf<-list(d,m,tw,Mx)
return(results.embed_udf)
}
