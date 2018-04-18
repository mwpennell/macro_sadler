library(castor)
library(diversitree)
source('R/ageGenerator.R')

ml_reg_stats_age_subtrees=function(phylo,ages,tolerance,N){
  lm.intercept<-c()
  lm.slope<-c()
 
 for (e in 1:N){
   new_age_generator<-ageGenerator(ages=ages,tree=phylo,tolerance=tolerance, fixed = T)
   new_age_generator$newGenera<-new_age_generator$newGenera[!sapply(new_age_generator$newGenera, anyNA)] 
   clade.age<-c()
   lambda<-c()
  
   for (i in 1:length(new_age_generator$newGenera)){
     clade.tree<-get_subtree_with_tips(phylo,new_age_generator$newGenera[[i]])$subtree
     if (Ntip(clade.tree)<=5) {next(i+1)} else 
     {clade.age[i]<-max(branching.times(clade.tree))
     l <- make.bd(clade.tree)
     f <- find.mle(l, x.init=c(.5,0.1))
     lambda[i]<- f$par["lambda"]}
   }  
  
   lm1<-lm(log(lambda)~log(clade.age))
   lm.intercept[[e]]<-lm1$coefficients[[1]]
   lm.slope[[e]]<-lm1$coefficients[[2]]
 }
  
  return(data.frame(lm.intercept, lm.slope))
}
