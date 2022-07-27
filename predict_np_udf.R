# A user-defined nonlinear prediction function
# After Huffaker, Bittelli & Rosa (2017) Nonlinear Time Series Analysis with R, 
# Oxford University Press, pp. 154-174.

predict_np_udf <- function(Mx,frac.learn){
  #Step 2a: Partition Mx into learning and test sets
  frac.learn<-0.5 #fraction in learning set
  learn.rows<-round(frac.learn*nrow(Mx))
  learn.em.0<-Mx[1:learn.rows,] #initial learning set
  test.em<-Mx[(learn.rows+1):nrow(Mx),] #initial test set
  #Step 2b: Prediction
  hold.test<-matrix(0,(nrow(test.em)),1)
  hold.pred<-matrix(0,(nrow(test.em)),1)
  for(i in 1:(nrow(test.em))) {
    #print("i");print(i)
    learn.em<-Mx[1:((nrow(learn.em.0)+i)-1),]
    #print("learn.em");print(learn.em)
    #Step 2b(1): Calculate nearest neighbours to last row in learning set (learn.em.0)
    #Distance between points on the attractor
    ref.point<-nrow(learn.em) #index of reference point
    library(fields)
    dist<-rdist(learn.em) #distance matrix
    dist.ref<-dist[,ref.point] #distances from reference point to other points
    sc<-max(dist.ref) #find maximum distance from reference point
    dist.ref.sc<-dist.ref/sc #scale distances from reference point to max distance
    #Order distances from closest to farthest from reference point
    o<-order(dist.ref.sc)
    #Remove distance of reference point to itself
    remove<-which(o==nrow(learn.em));o1<-o[-remove]
    #Indicies of m+1 smallest distances from reference point (m=embedding dimension)
    o2<-o1[1:(m+1)];o2<-na.omit(o2)
    dist.ordered<-dist.ref.sc[o2] #ordered distances
    #Increment neighbouring indices by 1 period for use in prediction algorithm
    o.pred<-o2+1
    pred.ngh<-learn.em[o.pred,] #Nearest neighbours
    #Step 2b(2): Compute prediction as average of neighbouring points weighted by
    # distances from reference point
    #Compute weights (Sughihara et al., 2012)
    u.denom<-dist.ordered[1] #distance from reference point to nearest neighbour
    hold<-matrix(0,(m+1),1)
    for(k in 1:(m+1)){ #summands in w.denom
      u.vector<-exp(-dist.ordered[k]/u.denom)
      hold[k,1]<-u.vector
    } #end loop k
    u.vector<-hold[1:(m+1)]
    w.denom<-sum(u.vector)
    w.vector<-u.vector/w.denom
    #Prediction of next point on attractor (row in test.em)
    pred.point<-w.vector%*%pred.ngh
    #Prediction of time series observation is first (unlagged) element of
    #pred.point
    pred.ts<-pred.point[1]
    #print("pred.ts");print(pred.ts)
    #Step 2b(3): Test point on attractor (row in test.em)
    test.point<-test.em[i,]
    #print("test.point");print(test.point)
    #Time series observation to be validated is first (unlagged) element of test.point
    test.ts<-test.point[1]
    hold.test[i,]<-test.ts #1st element is original data point
    hold.pred[i,]<-pred.ts
  } #end i loop through Mx
  #Step 3: Calculate modified Nash-Sutcliffe model efficiency (nse)
  num<-sum((hold.test-hold.pred)^2)
  den<-sum((hold.test-mean(hold.test))^2)
  nse<-1-(num/den)
  #print("nse");print(hold.mnse)
  #results<-cbind(learn.0,hold.test,hold.pred)
  #print("learn,test,pred");print(results)
  results<-list(nse,hold.test,hold.pred)
  return(results)
} #end function
results.np<-predict_np_udf(Mx,frac.learn=0.5)
nse<-results.np[[1]]
test<-results.np[[2]]
prediction<-results.np[[3]]
