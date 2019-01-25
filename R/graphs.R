library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(ggthemes)

#### Phylogenies ####
summary_tree_results<-read.delim("output/summary_tree_results.csv", sep = ",", dec=".")
express<-expression(paste("Ln mean ", lambda, " (species ", Myr^-1,")"),sep=" ")
e.express<-expression(paste("Ln mean ", mu, " (species ", Myr^-1,")"),sep=" ")

# Descriptive co-relations
a1<-ggplot(summary_tree_results, aes(x=log(tree.max.age),y=log(n.clade))) + geom_point() + 
  theme_tufte(base_family = "Helvetica") + geom_rangeframe() + theme(axis.title = element_text(size=15)) + 
  labs(x="Log Clade age (Myr)", y= "Log Clade richness") + 
  scale_x_continuous(breaks= c(2:6), labels=as.character(round(exp(c(2:6)),0))) +
  scale_y_continuous(breaks= seq(4,10,2), labels=as.character(round(exp(seq(4,10,2)),0))) +
  geom_smooth(method=lm, se=T,alpha=.2) 

a2<-ggplot(summary_tree_results, aes(x=log(tree.max.age),y=(gamma.stat))) + geom_point() + 
  theme_tufte(base_family = "Helvetica") + geom_rangeframe() +  theme(axis.title = element_text(size=15)) + 
  labs(x="Log Clade age (Myr)", y= "Gamma statistic") + 
  scale_x_continuous(breaks= c(2:6), labels=as.character(round(exp(c(2:6)),0))) +
  geom_smooth(method=lm, se=T,alpha=.2) 

ggarrange(a1, a2,  
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

g1<-ggplot(summary_tree_results, aes(x=(gamma.stat),y=log(mean.clade.lambda))) + geom_point() + 
  theme_tufte(base_family = "Helvetica") + geom_rangeframe() + theme(axis.title = element_text(size=15)) +
  labs(x= "Gamma statistic", y=express) + 
  geom_smooth(method=lm, se=T,alpha=.2)

g2<-ggplot(summary_tree_results, aes(x=(gamma.stat),y=log(mean.clade.mu))) + geom_point() +
  theme_tufte(base_family = "Helvetica") + geom_rangeframe() + theme(axis.title = element_text(size=15)) +
  labs(x= "Gamma statistic", y=e.express) + 
  geom_smooth(method=lm, se=T,alpha=.2)

ggarrange(g1, g2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

g3<-ggplot(summary_tree_results, aes(x=(gamma.stat),y=log(net.diver))) + geom_point() + 
  theme_tufte(base_family = "Helvetica") + geom_rangeframe() + theme(axis.title = element_text(size=15)) +
  labs(x= "Gamma statistic", y=expression(paste("Ln Net diversification (species ", Myr^-1,")"),sep=" ")) +  
  geom_smooth(method=lm, se=T,alpha=.2)

ggarrange(g1, g2, g3,  
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

# Fig. S1 Number of shifts
ss<-filter(summary_tree_results,summary_tree_results$best.n.shifts>=1)
hh1<-ggplot(ss, aes(x=best.n.shifts)) + geom_histogram(color="darkblue", fill="white") + 
  labs(title="", x="Number of Shifts", y = "") + theme_tufte(base_family = "Helvetica", base_size = 17) +
  geom_rangeframe(data=data.frame(x=c(0, 20), y=c(0, 80)), aes(x, y)) 

hh2<-ggplot(ss, aes(x=best.n.shifts/tree.max.age)) + geom_histogram(color="darkblue", fill="white") + 
  labs(title="",y = "", x=expression(paste("Number of Shifts per ", Myr^-1),sep=" "), y = "") + 
  theme_tufte(base_family = "Helvetica", base_size = 17) + 
  geom_rangeframe(data=data.frame(x=c(0, .4), y=c(0, 50)), aes(x, y)) 

ggarrange(hh1, hh2,   
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

round(1/mean(ss$best.n.shifts/ss$tree.max.age),1)
round(var(ss$best.n.shifts/ss$tree.max.age),3)

# Fig. 2A-B Time-dependency
ph1<-ggplot(summary_tree_results, aes(x=(tree.max.age),y=(mean.clade.lambda))) + 
  geom_point(show.legend = F, colour = "#0E233E", size=2) + theme_cowplot() + 
  labs(x="Clade age (Myr)", y=expression(paste("Speciation rate (species ", Myr^-1,")"),sep=" ")) + 
  theme(axis.title = element_text(size=15)) + 
  geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a=1,b=2)), se=F, na.rm = T, colour="#EA3770", size=1.5) + theme_tufte(base_family = "Helvetica") + 
  theme(axis.title = element_text(size=15)) +
  geom_rangeframe(data=data.frame(x = c(0, 350), y = c(0, 1.5)), aes(x, y))

ph2<-ggplot(summary_tree_results, aes(x=(tree.max.age),y=(mean.clade.mu))) + 
  geom_point(show.legend = F, colour = "#0E233E", size=2) + theme_cowplot() +
  labs(x="Clade age (Myr)", y=expression(paste("Extinction rate (species ", Myr^-1,")"),sep=" ")) +
  theme(axis.title = element_text(size=15)) +
  geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a=1,b=2)), se=F, na.rm = T, colour="#EA3770", size=1.5) + theme_tufte(base_family = "Helvetica") + 
  theme(axis.title = element_text(size=15)) +
  geom_rangeframe(data = data.frame(x = c(0, 350), y = c(0, .8)), aes(x, y))

l.inset<-ggplot(summary_tree_results, aes(x = log(tree.max.age), y = log(mean.clade.lambda))) + 
  geom_point(colour = "#0E233E") + theme_cowplot() + labs(x = "", y = "") + 
  scale_x_continuous(breaks = c(1,3.5,6), labels = as.character(round(exp(c(1,3.5,6)),0))) +
  scale_y_continuous(breaks = c(-4,-1.5,0.5), labels = as.character(round(exp(c(-4,-1.5,0.5)),2))) +
  # labs(x="Ln clade age (Myr)", y=express) +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20)) +
  geom_smooth(method=lm, se=T, alpha=.1, colour = "#EA3770") + theme_tufte(base_family = "Helvetica") +
  theme(axis.title = element_text(size = 20)) +
  geom_rangeframe(data=data.frame(x = c(1, 6), y = c(-4, 0.5)), aes(x, y))

m.inset<-ggplot(summary_tree_results, aes(x=log(tree.max.age),y=log(mean.clade.mu))) + 
  geom_point(colour = "#0E233E") + theme_cowplot() + labs(x = "", y = "") + 
  scale_x_continuous(breaks = c(1,3.5,6), labels = as.character(round(exp(c(1,3.5,6)),0))) +
  scale_y_continuous(breaks = c(-6,-3,0), labels = as.character(round(exp(c(-6,-3,0)),2))) +
  labs(x="", y="") +
  #labs(x="Ln clade age (Myr)", y=e.express) + 
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20)) +
  geom_smooth(method=lm, se=T,alpha=.1,colour="#EA3770") + theme_tufte(base_family = "Helvetica") +
  theme(axis.title = element_text(size = 20)) +
  geom_rangeframe(data=data.frame(x = c(1, 6), y = c(-6, 0)), aes(x, y))

gg1<-ggdraw() + draw_plot(ph1, x = 0, y = 0) + draw_plot(l.inset, x = 0.55, y = 0.6, width = 0.4, height = 0.4)
gg2<-ggdraw() + draw_plot(ph2, x = 0, y = 0) + draw_plot(m.inset, x = 0.55, y = 0.6, width = 0.4, height = 0.4)

ggarrange(gg1, gg2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

# Fig. S6 Null trees (push of the past)
empirical_bd_slope<-readRDS("null_trees/empirical/output/tree_bd_slope.rds")
m50_bd_slope<-readRDS("null_trees/mu50/output/tree_bd_slope.rds")
m75_bd_slope<-readRDS("null_trees/mu75/output/tree_bd_slope.rds")
null_bd_slope<-data.frame(empirical_bd_slope,m50_bd_slope,m75_bd_slope)

h1<-ggplot(null_bd_slope, aes(x=empirical_bd_slope)) + geom_histogram(color="darkblue", fill="white") + 
  geom_vline(aes(xintercept=-0.436),color="red", linetype="dashed", size=1) + 
  labs(title="Empirical parameters", x="Slope", y = "") + 
  theme_tufte(base_family = "Helvetica", base_size = 14) + 
  geom_rangeframe(data=data.frame(x=c(-0.436, .1), y=c(0, 150)), aes(x, y)) 

h2<-ggplot(null_bd_slope, aes(x=m50_bd_slope)) + geom_histogram(color="darkblue", fill="white") + 
  geom_vline(aes(xintercept=-0.436),color="red", linetype="dashed", size=1) + 
  labs(title="mu = 0.5 x Lambda",x="Slope", y = "") +
  theme_tufte(base_family = "Helvetica", base_size = 14) + 
  geom_rangeframe(data=data.frame(x=c(-0.436, .1), y=c(0, 155)), aes(x, y)) 

h3<-ggplot(null_bd_slope, aes(x=m75_bd_slope)) + geom_histogram(color="darkblue", fill="white") + 
  geom_vline(aes(xintercept=-0.436),color="red", linetype="dashed", size=1) + 
  labs(title="mu = 0.75 x Lambda",x="Slope", y = "") + 
  theme_tufte(base_family = "Helvetica", base_size = 14) + 
  geom_rangeframe(data=data.frame(x=c(-0.436, .1), y=c(0, 130)), aes(x, y)) 

ggarrange(h1, h2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

#### Fossils ####
long_paleo_results<-(gather(summary_paleo_results, key="Rate", value = "Estimate",mean.clade.mu:mean.clade.origination))
long_paleo_results$Rate<-recode(long_paleo_results$Rate, `mean.clade.mu`="Extinction", `mean.clade.origination`="Origination")

p1<-ggplot(summary_paleo_results, aes(x=(Duration), y=(mean.clade.origination), na.rm = TRUE)) + 
  geom_point(show.legend = F, colour = "#0E233E", size = 2) + theme_cowplot() + 
  labs(x="Duration (Myr)", y=expression(paste("Origination rate (genera ", Myr^-1,")"),sep=" ")) + 
  theme(axis.title = element_text(size = 15)) + 
  geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a=1,b=1)), se=F, colour="#EA3770", size=1.5) +
  theme_tufte(base_family = "Helvetica") + 
  theme(axis.title = element_text(size = 15)) +
  geom_rangeframe(data=data.frame(x = c(0, 500), y = c(0, .4)), aes(x, y))

p2<-ggplot(summary_paleo_results, aes(x = (Duration), y = (mean.clade.mu))) + 
  geom_point(show.legend = F, colour = "#0E233E", size=2) + theme_cowplot() +
  labs(x="Duration (Myr)", y=expression(paste("Extinction rate (genera ", Myr^-1,")"),sep=" ")) +
  theme(axis.title = element_text(size=15)) +
  geom_smooth(method = 'nls', formula = y~a*x^b, method.args = list(start = c(a=1,b=1)), se=F, na.rm = T, colour="#EA3770", size=1.5) +
  theme_tufte(base_family = "Helvetica") + 
  theme(axis.title = element_text(size = 15)) +
  geom_rangeframe(data = data.frame(x = c(0, 500), y = c(0, .8)), aes(x, y))

o.inset<-ggplot(summary_paleo_results, aes(x=log(Duration),y=log(mean.clade.origination))) + 
  geom_point(colour = "#0E233E") + theme_cowplot() + labs(x = "", y = "") + 
  scale_x_continuous(breaks = c(1,3.5,6), labels = as.character(round(exp(c(1,3.5,6)),0))) +
  scale_y_continuous(breaks = c(-4,-2.3, -.8), labels = as.character(round(exp(c(-4,-2.3,-.8)),2))) +
  # labs(x="Ln clade age (Myr)", y=express) +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20)) +
  geom_smooth(method=lm, se=T, alpha=.1, colour = "#EA3770") + theme_tufte(base_family = "Helvetica") +
  theme(axis.title = element_text(size = 20)) +
  geom_rangeframe(data=data.frame(x = c(1, 6), y = c(-4, -.8)), aes(x, y))

mf.inset<-ggplot(summary_paleo_results, aes(x=log(Duration),y=log(mean.clade.mu))) + 
  geom_point(colour = "#0E233E") + theme_cowplot() + labs(x = "", y = "") + 
  scale_x_continuous(breaks = c(1,3.5,6), labels = as.character(round(exp(c(1,3.5,6)),0))) +
  scale_y_continuous(breaks = c(-4,-1.5,0), labels = as.character(round(exp(c(-4,-1.5,0)),2))) +
  # labs(x="Ln clade age (Myr)", y=express) +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20)) +
  geom_smooth(method=lm, se=T, alpha=.1, colour = "#EA3770") + theme_tufte(base_family = "Helvetica") +
  theme(axis.title = element_text(size = 20)) +
  geom_rangeframe(data=data.frame(x = c(1, 6), y = c(-4, 0)), aes(x, y))

gg3<-ggdraw()+ draw_plot(p1, x = 0, y = 0) + draw_plot(o.inset, x = 0.55, y = 0.6, width = 0.4, height = 0.4)
gg4<-ggdraw()+ draw_plot(p2, x = 0, y = 0) + draw_plot(mf.inset, x = 0.55, y = 0.6, width = 0.4, height = 0.4)

#Fig. 2C-D Time-dependency
ggarrange(gg3, gg4, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

# Temporal trends
trend<-readRDS("output/Trends_pc.rds")
round(colMeans(trend[,2:5], na.rm = T),3)
apply(trend[,2:5],2, min, na.rm=T)

th1<-ggplot(trend, aes(x=Origination)) + geom_histogram(color="darkblue", fill="white") + 
  #geom_vline(aes(xintercept=-0.436),color="red", linetype="dashed", size=1) + 
  labs(title="", x="Spearman's rho", y = "") + 
  theme_tufte(base_family = "Helvetica", base_size = 14) + 
  geom_rangeframe(data=data.frame(x=c(-1, 1), y=c(0, 10)), aes(x, y)) 

th2<-ggplot(trend, aes(x=Extinction)) + geom_histogram(color="darkblue", fill="white") + 
  #geom_vline(aes(xintercept=-0.436),color="red", linetype="dashed", size=1) + 
  labs(title="", x="Spearman's rho", y = "") + 
  theme_tufte(base_family = "Helvetica", base_size = 14) + 
  geom_rangeframe(data=data.frame(x=c(-1, 1), y=c(0, 10)), aes(x, y)) 

th3<-ggplot(trend, aes(x=p.value.ori)) + geom_histogram(color="darkblue", fill="white") + 
  #geom_vline(aes(xintercept=-0.436),color="red", linetype="dashed", size=1) + 
  labs(title="", x="Spearman's rho", y = "") + theme_tufte(base_family = "Helvetica") + 
  geom_rangeframe(data=data.frame(x=c(0, 1), y=c(0, 10)), aes(x, y)) 

th4<-ggplot(trend, aes(x=p.value.ext)) + geom_histogram(color="darkblue", fill="white") + 
  #geom_vline(aes(xintercept=-0.436),color="red", linetype="dashed", size=1) + 
  labs(title="", x="Spearman's rho", y = "") + theme_tufte(base_family = "Helvetica") + 
  geom_rangeframe(data=data.frame(x=c(0, 1), y=c(0, 10)), aes(x, y)) 

ggarrange(th1, th2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

ggarrange(th3, th4, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

#### ####
#Fig. 1 Speciation - Origination rates histograms
express<-expression(paste("Mean ", lambda, " (species ", Myr^-1,")"),sep=" ")
express1<-expression(paste("Mean origination (genera ", Myr^-1,")"),sep=" ")

h1<-ggplot(summary_tree_results, aes(x=mean.clade.lambda)) +
  geom_histogram(colour="#0E233E", fill="white") + 
  labs(title="",x=express, y = "") + theme_tufte(base_family = "Helvetica", base_size = 24) +
  geom_rangeframe() + theme(axis.title = element_text(size = 24)) +
  geom_rangeframe(data=data.frame(x=c(0, 1.6), y=c(0, 30)), aes(x, y))

h2<-ggplot(summary_paleo_results, aes(x=mean.clade.origination)) +
  geom_histogram(colour="#0E233E", fill="white") + 
  labs(title="",x=express1, y = "") + theme_tufte(base_family = "Helvetica", base_size = 24) + 
  geom_rangeframe() + theme(axis.title = element_text(size = 24)) + 
  geom_rangeframe(data=data.frame(x=c(0, 0.5), y=c(0, 20)), aes(x, y))

ggarrange(h1, h2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

# Fig 2. Time-dependency of macroevolutionary
ggarrange(gg1, gg2, gg3, gg4, 
          labels = c("A", "B", "C", "D"),
          font.label = list(size = 16),
          ncol = 2, nrow = 2)

#Extinction rates histograms
e.express<-expression(paste("Mean ", mu, " (species ", Myr^-1,")"),sep=" ")
e.express1<-expression(paste("Mean extinction (genera ", Myr^-1,")"),sep=" ")

hh1<-ggplot(summary_tree_results, aes(x=mean.clade.mu)) + 
  geom_histogram(colour="#0E233E", fill="white") + 
  labs(title="",x=e.express, y = "") + theme_tufte(base_family = "Helvetica",  base_size = 14) + 
  geom_rangeframe() + 
  theme(axis.title = element_text(size=15)) + 
  geom_rangeframe(data=data.frame(x=c(0, 1), y=c(0, 30)), aes(x, y))

hh2<-ggplot(summary_paleo_results, aes(x=mean.clade.mu)) + 
  geom_histogram(colour="#0E233E", fill="white") + 
  labs(title="",x=e.express1, y = "") + theme_tufte(base_family = "Helvetica",  base_size = 14) +
  geom_rangeframe() + 
  theme(axis.title = element_text(size=15)) + 
  geom_rangeframe(data=data.frame(x=c(0, 1), y=c(0, 30)), aes(x, y))

ggarrange(hh1, hh2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
