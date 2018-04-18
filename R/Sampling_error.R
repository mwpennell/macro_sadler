library(diversitree)
library(rncl)

# Load Jetz et al. 2012 bird family trees
j.phylo.path<-("/Users/luisfranciscohenaodiaz/Documents/UBC/Projects/spp_rates/phylosOrig/Jetz et al. /MCC_trees/data/")

d<-dir(j.phylo.path)
t<-grep("tree",d)
tt<-d[t]
j.phylos<-list()
for(i in 1:length(tt)){
  tree<-read_newick_phylo(paste0(j.phylo.path,tt[i]))
  if(Ntip(tree)>=10)
    {j.phylos[[i]]<-tree} else {"NA"}
}
j.phylos<-j.phylos[!sapply(j.phylos, is.null)] 
head(j.phylos)

# bird family tree ages
j.ages<-c()
for (i in 1:length(j.phylos)){
  j.ages[i]<-max(branching.times(j.phylos[[i]]))
}

source('R/ML_birth-death.R', chdir = TRUE)
j.ml<-ml_bd.stats_trees(c(j.phylos[1:27],j.phylos[29:38],j.phylos[40:46]))

# Load Zanne et al. plant order trees
z.phylo.path<-("/Users/luisfranciscohenaodiaz/Documents/UBC/Projects/spp_rates/phylosOrig/Zanne et al./Zanne Pruned/")
d<-dir(z.phylo.path)
t<-grep("tree",d)
tt<-d[t]

z.phylos<-list()
for(i in 1:length(tt)){
  tree<-read_newick_phylo(paste0(z.phylo.path,tt[i]))
  if(Ntip(tree)>=10)
  {z.phylos[[i]]<-tree} else {"NA"}
}
z.phylos<-z.phylos[!sapply(z.phylos, is.null)] 
head(z.phylos)

# plant order tree ages
z.ages<-c()
for (i in 1:length(z.phylos)){
  z.ages[i]<-max(branching.times(z.phylos[[i]]))
}

z.ml<-ml_bd.stats_trees(z.phylos)

###
j.phylo<-read_newick_phylo("/Users/luisfranciscohenaodiaz/Documents/UBC/Projects/spp_rates/phylosOrig/tree_birdtree_.txt")

z.phylo<-read_newick_phylo("/Users/luisfranciscohenaodiaz/Documents/UBC/Projects/spp_rates/phylosOrig/tree_Z2014_.tre")
z.phylo<-castor:::extend_tree_to_height(z.phylo)$tree
z.phylo$node.label<-NULL

source('R/ml_reg_stats_age_subtrees.R', chdir = TRUE)

j.ml.chop<-ml_reg_stats_age_subtrees(j.phylo,j.ages,tolerance=10,N=100)
z.ml.chop<-ml_reg_stats_age_subtrees(z.phylo,z.ages,tolerance=10,N=100)

par(mfrow=c(1,2))
hist(j.ml.chop$lm.slope, main="Bird families",xlab = "Slope")
abline(v=coefficients(lm(j.ml$l.lambda~j.ml$clade.age))[2], col="red")
coefficients(lm(j.ml$l.lambda~j.ml$clade.age))[2]

hist(z.ml.chop$lm.slope, main="Plant orders",xlab = "Slope")
abline(v=coefficients(lm(z.ml$l.lambda~z.ml$clade.age))[2], col="red")
coefficients(lm(z.ml$l.lambda~z.ml$clade.age))[2]

##
tree_metrics=function(tt){
  phylogs<-c()
  ntips<-c()
  tree.min.age<-c()
  tree.max.age<-c()
  gamma.stat<-c()
  trees.metrics<-c()
  for (i in 1:length(tt)){
    phylogs<-tt[[i]]
    ntips[i]<-Ntip(phylogs)
    tree.min.age[i]<-min(branching.times(phylogs))
    tree.max.age[i]<-max(branching.times(phylogs))
    gamma.stat[i]<-gammaStat(phylogs)
    trees.metrics<-data.frame(ntips,tree.min.age,tree.max.age,gamma.stat)
  }  
  return(trees.metrics)
}

j.t.metrics<-tree_metrics(j.phylos)
z.t.metrics<-tree_metrics(z.phylos)

par(mfrow=c(1,2))
hist(j.t.metrics$tree.max.age,main="Bird families",xlab="Ages")
hist(z.t.metrics$tree.max.age,main="Plant orders",xlab="Ages")
