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

j.ml.chop<-ml_reg_stats_age_subtrees(j.phylo,j.ages,tolerance=10,N=1000)
z.ml.chop<-ml_reg_stats_age_subtrees(z.phylo,z.ages,tolerance=10,N=1000)
t.ml.chop<-ml_reg_stats_age_subtrees(t.phylo,t.ages,tolerance=10,N=1000)

## Figure 4. Random clade selection effect
h1<-ggplot(j.ml.chop, aes(x=lm.slope)) + geom_histogram(colour= "#0E233E", fill="white") + 
  geom_vline(aes(xintercept=-0.436), colour="#EA3770", linetype="dashed", size=1) + 
  labs(title="", x="Slope", y = "") +
  theme_tufte(base_family = "Helvetica", base_size = 15) + 
  geom_rangeframe(data=data.frame(x=c(-3, 1), y=c(0, 300)), aes(x, y)) 

h2<-ggplot(z.ml.chop, aes(x=lm.slope)) + geom_histogram(color="#0E233E", fill="white") + 
  geom_vline(aes(xintercept=-0.436),color="#EA3770", linetype="dashed", size=1) + 
  labs(title="",x="Slope", y = "") + theme_tufte(base_family = "Helvetica", base_size = 15) +
  geom_rangeframe(data=data.frame(x=c(-3, 1), y=c(0, 150)), aes(x, y)) 

h3<-ggplot(t.ml.chop, aes(x=lm.slope)) + geom_histogram(color="#0E233E", fill="white") + 
  geom_vline(aes(xintercept=-0.436),color="#EA3770", linetype="dashed", size=1) + 
  labs(title="",x="Slope", y = "") + theme_tufte(base_family = "Helvetica", base_size = 15) +
  geom_rangeframe(data=data.frame(x=c(-2, 1), y=c(0, 150)), aes(x, y)) 

ggarrange(h1, h2, h3, 
          labels = c("A", "B", "C"),
          font.label = list(font.label = 20),
          ncol = 3, nrow = 1)
 