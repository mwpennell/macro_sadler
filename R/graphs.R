library(ggplot2)
library(RColorBrewer)

summary.tree.results<-read.delim("output/summary_tree_results.csv", sep = ",", dec=".")

express<-expression(paste("Ln mean ", lambda, " (species ", Myr^-1,")"),sep=" ")
g<-ggplot(summary.tree.results, aes(x=log(tree.max.age),y=log(mean.clade.lambda)))+ 
  geom_point(aes(size=n.clade,colour=gamma.stat))+theme_bw() +scale_size_continuous(name = "Clade\nrichness") 

g+labs(x="Ln clade age (Myr)", y=express) + theme(axis.title = element_text(size=15)) + 
  geom_smooth(method=lm, se=T,alpha=.1)  + scale_colour_gradientn(name="Gamma\nstatistic",colours = rev(brewer.pal(5,"Spectral"))) + guides(size = guide_legend(order=1))

summary.paleo.results<-read.delim("output/summary_paleo_results.csv", sep = ",", dec=".")

g<-ggplot(summary.paleo.results,aes(x=log(Duration),y=log(mean.clade.origination)))+ geom_point(aes(size=Ngen,colour=freqRat))+theme_bw()+scale_size_continuous(name = "Clade\nrichness") 

g+labs(x="Ln duration (Myr)", y="Estimate")+ scale_colour_gradientn(name="Sampling\ncompletedness",colours = rev(brewer.pal(5,"Spectral"))) + guides(size = guide_legend(order=1))
