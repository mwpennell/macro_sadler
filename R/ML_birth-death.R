library(diversitree)
library(visreg)

ml_bd.stats_trees=function(sim.trees){
  clade.size<-c()
  clade.age<-c()
  lambda<-c()
  mu<-c()
  for(i in 1:length(sim.trees)){
      clade.tree<-sim.trees[[i]]
      clade.size[i]<-length(clade.tree$tip.label)
      clade.age[i]<-max(branching.times(clade.tree))
      l <- make.bd(clade.tree)
      f <- find.mle(l, x.init=c(0.5,0.1))
      lambda[i]<- f$par["lambda"]
      mu[i]<- f$par["mu"]
      sum.tree<-data.frame(clade.size,clade.age,lambda,mu,l.clade.size=log(clade.size), l.clade.age=log(clade.age), l.lambda=log(lambda), l.mu=log(mu))
  }
  return(sum.tree)
}

ml_bd.constrain.stats_trees=function(sim.trees){
  clade.size<-c()
  clade.age<-c()
  lambda<-c()
  for(i in 1:length(sim.trees)){
    clade.tree<-sim.trees[[i]]
    clade.size[i]<-length(clade.tree$tip.label)
    clade.age[i]<-max(branching.times(clade.tree))
    l <- make.bd(clade.tree)
    l2 <-constrain(l, mu~0)
    f <- find.mle(l2, x.init=c(1))
    lambda[i]<- f$par["lambda"]
    sum.tree<-data.frame(clade.size,clade.age,lambda,l.clade.size=log(clade.size), l.clade.age=log(clade.age), l.lambda=log(lambda))
  }
  return(sum.tree)
}

