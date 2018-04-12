## compute support for time-dependent models
library(diversitree)
library(geiger)
compare_bd_t <- function(tree, sampling.f){
  lik <- make.bd(tree, sampling.f=sampling.f)
  fit <- find.mle(lik, x.init=c(0.5, 0.1))
  lik2 <- make.bd.t(tree, sampling.f=sampling.f,
                    functions=c("linear.t", "constant.t"))
  fit2 <- find.mle(lik2, x.init=c(0.5, 0.05, 0.1))
  aic_bd <- AIC(fit)
  aic_bdt <- AIC(fit2)
  aic_w <- aicw(c(aic_bd, aic_bdt))
  list(slope_est = as.numeric(fit2$par["lambda.m"]), aicw_t_model = aic_w[2, "w"])
}

phy.path<-("data/phylo_")

d<-dir("data")
t<-grep("phylo_",d)
tt<-d[t]
emp.phylogs<-list()
for(i in 1:length(tt)){
  tree <- read.tree(paste0("data/",tt[[i]]))
  if(Ntip(tree)>=10)
  {emp.phylogs[[i]]<-tree} else {"NA"}
  }
emp.phylogs<-emp.phylogs[!sapply(emp.phylogs, is.null)] 

summary_tree_results<-read.csv("output/summary_tree_results.csv",header=T,row.names = 1)
summary_tree_results<-subset(summary_tree_results, summary_tree_results$ntips>=10)

slope<-c()
aicw_t_model<-c()
res<-c()
for(i in 1:length(emp.phylogs)){
  slope[i]<-compare_bd_t(emp.phylogs[[i]],summary_tree_results$prop.samp[i])$slope_est
  aicw_t_model[i]<-compare_bd_t(emp.phylogs[[i]],summary_tree_results$prop.samp[i])$aicw_t_model
  res<-data.frame(slope,aicw_t_model)
}
res
summary(res)

write.table(res,file = "output/time_dependence.csv", sep=",", dec = ".")
