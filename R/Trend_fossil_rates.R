

library(paleotree)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(ggplot2)
library(corrplot)

results = readRDS("data/Trends_pc.rds")
head(results)

res2gg = results %>%
  select(-Diversification, -p.value.div, -Duration) %>%
  gather(rate, value, -order, -starts_with("p.value."))


sm = ggplot(res2gg, aes(value)) + # , y = ..density..
  geom_histogram(color="darkblue", fill="white") +
  xlim(-1, 1) + 
  labs(title="",y = "Count", x = "Spearman's rho") + 
  theme_tufte(base_family = "Helvetica") + 
  geom_rangeframe(data=data.frame(x=c(-1, 1), y=c(0, 10)), aes(x, y)) +
  facet_grid(. ~ rate)
sm
ggsave(filename = "output/Trend_rates_percapita.png", plot = sm,
       width = 20, height = 14, units = "cm")

# Correlation as function of clade duration
trend = ggplot(results, aes(x = Duration, y = Diversification)) +
  geom_point(colour= "#0E233E") +
  theme_tufte(base_family = "Helvetica") +
  geom_rangeframe() +
  labs(title = "", x = "Clade age (Myr)", y = "Spearman's correlation (rho)") +
  geom_hline(aes(yintercept = 0), colour = "#EA3770", linetype = "dashed", size = 1) + 
  theme(axis.title = element_text(size=15)) + 
  #  geom_smooth(method = "lm", color = "black") +
  geom_rangeframe(data = data.frame(x = c(0, 500), y = c(-1,1)), aes(x, y))
trend
ggsave(filename = "output/Trend_rates_percapita_time.png", plot = trend)

model = lm(Diversification ~ Duration, data = results)
anova(model)
#

