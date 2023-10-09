library("tidyr")
library("dplyr")
library("ggplot2")
library("patchwork")

## -----------------------------------------
## Create the Plots
## -----------------------------------------

load("res.RData")

# Aggregate results from all three settings
results_all <-
  bind_rows(
    results_binom %>% as_tibble %>% mutate(setting = "Binomial PCA"),
    results_gammasmall %>% as_tibble %>% mutate(setting = "Gamma Clustering - Small"),
    results_gammalarge %>% as_tibble %>% mutate(setting = "Gamma Clustering - Large")
  ) %>%
  mutate(across(c("sim", "k", "value"), ~as.numeric(.)))

# Create the NLL as a function of k plot (Figure 3)
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

pvalNLL

# Create the NLL histogram (Figure 5)
phistNLL <- results_all %>% 
  filter(measure == "NLL") %>%
  group_by(setting, sim, method) %>%
  filter(value == min(value)) %>%
  ungroup %>%
  filter(method %in% c("Data Thinning (0.8)", "Multifold Thinning")) %>%
  ggplot(aes(x=k, fill=method)) +
  geom_vline(data=data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
                             val = c(10, 4, 10)), aes(xintercept=val)) +
  geom_bar(aes(y = (..count..)/nreps), position=position_dodge2(preserve="single")) + 
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


phistNLL

# Create the MSE as a function of k plot (Figure 7)
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

pvalMSE

# Create the MSE histogram (Figure 9)
phistMSE <- results_all %>% 
  filter(measure == "MSE") %>%
  group_by(setting, sim, method) %>%
  filter(value == min(value)) %>%
  ungroup %>%
  filter(method %in% c("Data Thinning (0.8)", "Multifold Thinning")) %>%
  ggplot(aes(x=k, fill=method)) +
  geom_vline(data=data.frame(setting = c("Binomial PCA", "Gamma Clustering - Small", "Gamma Clustering - Large"),
                             val = c(10, 4, 10)), aes(xintercept=val)) +
  geom_bar(aes(y = (..count..)/nreps), position=position_dodge2(preserve="single")) + 
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

phistMSE
