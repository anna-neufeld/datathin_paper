library(tidyverse)
library(patchwork)
setwd("~/datathin_paper/Figures_4_6_8_10/resAD/res")


#### READ IN SMALL GAMMA RES. 
file_names <- dir('.', pattern="dec16_gammasmall") #where you have your files

results_gammasmall <- read.csv(file_names[1], sep="", header=FALSE)
print(nrow(results_gammasmall))
names(results_gammasmall) <- c("eps", "rand", "lucyLoss", "sim", "method", "measure", "k", "value")
counter <- max(results_gammasmall$sim)
for (f in file_names[-1]) {
  temp <- read.csv(f, sep="", header=FALSE)
  print(nrow(temp))
  names(temp) <- c("eps", "rand", "lucyLoss", "sim", "method", "measure", "k", "value")
  temp$sim <- temp$sim+counter 
  counter <- max(temp$sim)
  results_gammasmall <- rbind(results_gammasmall, temp)
}

#### READ IN GAMMA LARGE
file_names <- dir(".", pattern="dec15_gammabig") #where you have your files

results_gammalarge <- read.csv(file_names[1], sep="", header=FALSE)
print(nrow(results_gammalarge))
names(results_gammalarge) <- c("eps", "rand", "lucyLoss", "sim", "method", "measure", "k", "value")
counter <- max(results_gammalarge$sim)
for (f in file_names[-1]) {
  temp <- read.csv(f, sep="", header=FALSE)
  print(nrow(temp))
  names(temp) <- c("eps", "rand", "lucyLoss", "sim", "method", "measure", "k", "value")
  temp$sim <- temp$sim+counter 
  counter <- max(temp$sim)
  results_gammalarge <- rbind(results_gammalarge, temp)
}


#### READ IN BINOM
file_names <- dir(".", pattern="dec15_bin") #where you have your files

results_binom <- read.csv(file_names[1], sep="", header=FALSE)
print(nrow(results_binom))
names(results_binom ) <- c("eps", "sim", "measure", "k", "value")
counter <- max(results_binom$sim)
for (f in file_names[-1]) {
  temp <- read.csv(f, sep="", header=FALSE)
  print(nrow(temp))
  names(temp) <- c("eps", "sim", "measure", "k", "value")
  temp$sim <- temp$sim+counter 
  counter <- max(temp$sim)
  results_binom <- rbind(results_binom, temp)
}



consRes_gammasmall_pre <- results_gammasmall %>% group_by(sim, eps, method, measure) %>%
  summarize(minK = which.min(value),
            bestRand = which.max(rand),
  )
consRes_gammasmall <- consRes_gammasmall_pre %>% group_by(eps, method, measure) %>% 
  summarize(propCorrect = mean(minK==4))


pgamma_small <- ggplot(data=consRes_gammasmall %>% filter(measure=="NLL"), aes(x=eps, y=propCorrect))+geom_point()+geom_line()+
 ggtitle("Proportion of times we select correct K")+
  theme_bw()+
  xlab(expression(epsilon))+ylab("Proportion of simulated datasets with correct K")+
  ylim(c(0,1)) +
  ggtitle("Gamma Clustering - Small")


consRes_gammalarge_pre <- results_gammalarge %>% group_by(sim, eps, method, measure) %>%
  summarize(minK = which.min(value),
            bestRand = which.max(rand),
  )


consRes_gammalarge <- consRes_gammalarge_pre %>% group_by(eps, method, measure) %>% 
  summarize(propCorrect = mean(minK==10), bestRand = mean(bestRand==10))


pgamma_large <- ggplot(data=consRes_gammalarge %>% filter(measure=="NLL"), aes(x=eps, y=propCorrect))+geom_point()+geom_line()+
  ggtitle("Proportion of times we select correct K")+
  theme_bw()+
  xlab(expression(epsilon))+ylab("Proportion of simulated datasets with correct K")+
  ylim(c(0,1)) +
  ggtitle("Gamma Clustering - Large")


consRes_binom_pre <- results_binom %>% group_by(sim, eps, measure) %>%
  summarize(minK = which.min(value))


consRes_binom <- consRes_binom_pre %>% group_by(eps,measure) %>% 
  summarize(propCorrect = mean(minK==10))


pbinom <- ggplot(data=consRes_binom %>% filter(measure=="NLL"), aes(x=eps, y=propCorrect))+geom_point()+geom_line()+
  theme_bw()+
  xlab(expression(epsilon^{(train)}))+ylab("Proportion of simulated datasets with correct K")+
  ylim(c(0,1)) +
  ggtitle("Binomial PCA")


pgamma_small2 <- ggplot(data=consRes_gammasmall %>% filter(measure=="MSE"), aes(x=eps, y=propCorrect))+geom_point()+geom_line()+
  ggtitle("Proportion of times we select correct K")+
  theme_bw()+
  xlab(expression(epsilon))+ylab("Proportion of simulated datasets with correct K")+
  ylim(c(0,1)) +
  ggtitle("Gamma Clustering - Small")

pgamma_large2 <- ggplot(data=consRes_gammalarge %>% filter(measure=="MSE"), aes(x=eps, y=propCorrect))+geom_point()+geom_line()+
  ggtitle("Proportion of times we select correct K")+
  theme_bw()+
  xlab(expression(epsilon))+ylab("Proportion of simulated datasets with correct K")+
  ylim(c(0,1)) +
  ggtitle("Gamma Clustering - Large")

pbinom2 <- ggplot(data=consRes_binom %>% filter(measure=="RealMSE"), aes(x=eps, y=propCorrect))+geom_point()+geom_line()+
  theme_bw()+
  xlab(expression(epsilon))+ylab("Proportion of simulated datasets with correct K")+
  ylim(c(0,1)) +
  ggtitle("Binomial PCA")




#### READ IN ALL METHODS RESULTS. 
file_names <- dir(".", pattern="dec16_allsims") #where you have your files
file_names <- file_names[grep("binomial", file_names)]

results_allbin <- read.csv(file_names[1], sep="", header=FALSE)
names(results_allbin) <- c("method", "sim", "measure", "k", "value")
print(nrow(results_allbin))
counter <- max(results_allbin$sim)
for (f in file_names[-1]) {
  temp <- read.csv(f, sep="", header=FALSE)
  print(nrow(temp))
  names(temp) <- c("method", "sim", "measure", "k", "value")
  temp$sim <- temp$sim+counter 
  counter <- max(temp$sim)
  results_allbin <- rbind(results_allbin, temp)
}

file_names <- dir(".", pattern="dec16_allsims") #where you have your files
file_names <- file_names[grep("gammasmall", file_names)]

results_allgs <- read.csv(file_names[1], sep="", header=FALSE)
names(results_allgs) <- c("method", "rand", "sim", "measure", "k", "value")
print(nrow(results_allgs))
counter <- max(results_allgs$sim)
for (f in file_names[-1]) {
  temp <- read.csv(f, sep="", header=FALSE)
  print(nrow(temp))
  names(temp) <- c("method", "rand", "sim", "measure", "k", "value")
  temp$sim <- temp$sim+counter 
  counter <- max(temp$sim)
  results_allgs <- rbind(results_allgs, temp)
}

file_names <- dir(".", pattern="dec16_allsims") #where you have your files
file_names <- file_names[grep("gammabig", file_names)]

results_allgb <- read.csv(file_names[1], sep="", header=FALSE)
names(results_allgb) <- c("method", "rand", "sim", "measure", "k", "value")
print(nrow(results_allgb))
counter <- max(results_allgb$sim)
for (f in file_names[-1]) {
  temp <- read.csv(f, sep="", header=FALSE)
  print(nrow(temp))
  names(temp) <- c("method", "rand", "sim", "measure", "k", "value")
  temp$sim <- temp$sim+counter 
  counter <- max(temp$sim)
  results_allgb <- rbind(results_allgb, temp)
}

results_all <-
  bind_rows(
    results_allbin %>% mutate(setting = "Binomial PCA", measure = ifelse(measure == "RealMSE", "MSE", measure)),
    results_allgs %>% select(-rand) %>% mutate(setting = "Gamma Clustering - Small"),
    results_allgb %>% select(-rand) %>% mutate(setting = "Gamma Clustering - Large")
  ) %>%
  mutate(
    method = case_when(
      method == "DataThin0.5" ~ "Data Thinning (0.5)",
      method == "DataThin0.8" ~ "Data Thinning (0.8)",
      method == "MultifoldThin" ~ "Multifold Thinning",
      TRUE ~ method
  ))

pvalNLL <- results_all %>%
  filter(measure == "NLL") %>%
  group_by(setting, method, k) %>%
  summarise(value = mean(value)) %>%
  mutate(value = (value-min(value))/(max(value)-min(value))) %>%
  ungroup %>%
  ggplot(aes(x=k, y=value, colour=method)) +
  geom_line() +
  geom_point(data=(results_all %>% 
                     filter(measure == "NLL") %>%
                     group_by(setting, method, k) %>% 
                     summarise(value = mean(value)) %>% 
                     mutate(value = (value-min(value))/(max(value)-min(value))) %>%
                     filter(value == min(value))
  ), 
  pch=1, size=3, stroke=2) +
  facet_wrap(~factor(setting, levels=c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large")), scales="free_x") +
  scale_colour_manual(name=NULL, 
                      values=c("Data Thinning (0.5)" = "#F8766D",
                               "Data Thinning (0.8)" = "#00BFC4",
                               "Multifold Thinning" = "#7CAE00",
                               "Naive" = "#C77CFF")) +
  geom_vline(data=data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
                             val = c(10, 4, 10)), aes(xintercept=val)) +
  xlab("K") + 
  ylab("") +
  theme(legend.position = "right") +
  theme_bw() +
  coord_cartesian(ylim=c(0,0.25))


pvalMSE <- results_all %>%
  filter(measure == "MSE") %>%
  group_by(setting, method, k) %>%
  summarise(value = mean(value)) %>%
  mutate(value = (value-min(value))/(max(value)-min(value))) %>%
  ungroup %>%
  ggplot(aes(x=k, y=value, colour=method)) +
  geom_line() +
  geom_point(data=(results_all %>% 
                     filter(measure == "MSE") %>%
                     group_by(setting, method, k) %>% 
                     summarise(value = mean(value)) %>% 
                     mutate(value = (value-min(value))/(max(value)-min(value))) %>%
                     filter(value == min(value))
  ), 
  pch=1, size=3, stroke=2) +
  facet_wrap(~factor(setting, levels=c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large")), scales="free_x") +
  scale_colour_manual(name=NULL, 
                      values=c("Data Thinning (0.5)" = "#F8766D",
                               "Data Thinning (0.8)" = "#00BFC4",
                               "Multifold Thinning" = "#7CAE00",
                               "Naive" = "#C77CFF")) +
  geom_vline(data=data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
                             val = c(10, 4, 10)), aes(xintercept=val)) +
  xlab("K") + 
  ylab("") +
  theme(legend.position = "right") +
  theme_bw() +
  coord_cartesian(ylim=c(0,0.25))


errdat <-  results_all %>% 
  filter(measure == "NLL") %>%
  group_by(setting, sim, method) %>%
  filter(value == min(value)) %>%
  ungroup %>%
  filter(method %in% c("Data Thinning (0.8)", "Multifold Thinning")) %>% 
  group_by(method, measure, setting, k) %>% 
  summarise(n=n()) %>% 
  mutate(p = n/sum(n),
         l = p - sqrt(p*(1-p)/2000),
         u = p + sqrt(p*(1-p)/2000))

phistNLL <- results_all %>% 
  filter(measure == "NLL") %>%
  group_by(setting, sim, method) %>%
  filter(value == min(value)) %>%
  ungroup %>%
  filter(method %in% c("Data Thinning (0.8)", "Multifold Thinning")) %>%
  ggplot(aes(x=k, fill=method)) +
  geom_vline(data=data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
                             val = c(10, 4, 10)), aes(xintercept=val)) +
  geom_bar(aes(y = (..count..)/2000), position=position_dodge2(preserve="single")) + 
  # geom_errorbar(data=errdat, aes(x=k, ymin=l, ymax=u), position=position_dodge2(preserve="single")) +
  facet_wrap(~factor(setting, levels=c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large")), scales="free_x") +
  scale_fill_manual(name=NULL, 
                    values=c("Data Thinning (0.8)" = "#00BFC4",
                             "Multifold Thinning" = "#7CAE00")) +
  geom_blank(data=expand_grid(
    data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
               val = c(6, 1, 6)),
    data.frame(method = c("Data Thinning (0.8)", "Multifold Thinning"))),
    aes(x=val)) +
  geom_blank(data=expand_grid(
    data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
               val = c(16, 11, 16)),
    data.frame(method = c("Data Thinning (0.8)", "Multifold Thinning"))),
    aes(x=val)) +
  ylim(c(0,1)) +
  xlab("K") + 
  ylab("") +
  theme(legend.position = "right") +
  theme_bw()



phistMSE <- results_all %>% 
  filter(measure == "MSE") %>%
  group_by(setting, sim, method) %>%
  filter(value == min(value)) %>%
  ungroup %>%
  filter(method %in% c("Data Thinning (0.8)", "Multifold Thinning")) %>%
  ggplot(aes(x=k, fill=method)) +
  geom_vline(data=data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
                             val = c(10, 4, 10)), aes(xintercept=val)) +
  geom_bar(aes(y = (..count..)/2000), position=position_dodge2(preserve="single")) + 
  facet_wrap(~factor(setting, levels=c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large")), scales="free_x") +
  scale_fill_manual(name=NULL, 
                    values=c("Data Thinning (0.8)" = "#00BFC4",
                             "Multifold Thinning" = "#7CAE00")) +
  geom_blank(data=expand_grid(
    data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
               val = c(6, 1, 6)),
    data.frame(method = c("Data Thinning (0.8)", "Multifold Thinning"))),
    aes(x=val)) +
  geom_blank(data=expand_grid(
    data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
               val = c(16, 11, 16)),
    data.frame(method = c("Data Thinning (0.8)", "Multifold Thinning"))),
    aes(x=val)) +
  ylim(c(0,1)) +
  xlab("K") + 
  ylab("") +
  theme(legend.position = "right") +
  theme_bw()



###### SAVE THE ACTUAL FINAL PLOTS
setwd("~/Dropbox/Generalized Count Splitting/JMLR-resubmit-sep-2023-v2/Figures"
      )
pvalNLL& 
  theme(axis.text = element_text(size=14), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), plot.title=element_text(size=16), legend.text=element_text(size=16),
        strip.text=element_text(size=16)) & ylab("Normalized NLL Loss")
ggsave("NLLcurves.png", width=12, height=4,  dpi=600) 

pvalMSE& 
  theme(axis.text = element_text(size=14), axis.title.x=element_text(size=20), plot.title=element_text(size=16), axis.title.y=element_text(size=20), legend.text=element_text(size=16),
        strip.text=element_text(size=16)) & ylab("Normalized MSE Loss")
ggsave("MSEcurves.png", width=12, height=4,  dpi=600) 


pbinom+pgamma_small+pgamma_large & 
  xlab(expression(epsilon^{(train)})) & theme(axis.text = element_text(size=14), axis.title.y=element_text(size=14), axis.title.x=element_text(size=20), plot.title=element_text(size=16))
ggsave("role_of_eps2.png", width=12, height=5,  dpi=600) 

pbinom2+pgamma_small2+pgamma_large2 & xlab(expression(epsilon^{(train)})) & theme(axis.text = element_text(size=14), axis.title.y=element_text(size=14), axis.title.x=element_text(size=20), plot.title=element_text(size=16))
ggsave("role_of_eps3.png", width=12, height=5,  dpi=600)

phistMSE + theme(axis.text = element_text(size=14), axis.title.x=element_text(size=20), axis.title.y=element_text(size=16), plot.title=element_text(size=16),
                 strip.text=element_text(size=16), legend.text=element_text(size=16))&ylab("Proportion of simulated datasets")
ggsave("MSEhist.png", width=12, height=4,  dpi=600)

phistNLL + 
  theme(axis.text = element_text(size=14), axis.title.x=element_text(size=20), axis.title.y=element_text(size=16), plot.title=element_text(size=16), 
        strip.text=element_text(size=16), legend.text=element_text(size=16)) & ylab("Proportion of simulated datasets")

ggsave("NLLhist.png", width=12, height=4,  dpi=600)




