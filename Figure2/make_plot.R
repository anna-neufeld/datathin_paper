library(tidyverse)
library(patchwork)

setwd("~/datathin_paper/Figure2")
load("fixedX_res-nov6-1out.RData")
load("randomX_res-nov6.RData")

res.FIXED$setting <-  "Non iid" 
res.RANDOM$setting <- "iid"
res <- rbind(res.FIXED, res.RANDOM)

prop.select <- function(u) {mean(!is.na(u))}

res.detect <- res %>% group_by(setting, method, beta, epsilon) %>% summarize_all(prop.select) %>% mutate(eps2 = paste("\u03B5 = ", epsilon))
res.power <- res %>% group_by(setting, method, beta, epsilon) %>% summarize_all(mean, na.rm=T) %>% mutate(eps2 = paste("\u03B5 = ", epsilon))

p1 <- ggplot(data=res.detect)+
  geom_line(lwd=1.3, aes(x=beta, y=X3, col=method, lty=eps2))+ggtitle("Detection")+ylab("Proportion of simulated datasets")
p2 <- ggplot(data=res.power, aes(x=beta, y=X3, col=method, lty=eps2))+geom_line(lwd=1.3)+ylab("Proportion of simulated datasets")+ggtitle("Power")
p1+p2 + plot_layout(guides="collect", ncol=3) &
  scale_linetype_manual(values=c(1,3,10)) & 
  theme_bw() & facet_grid(vars(setting)) &
  theme(axis.text=element_text(size=16), strip.text=element_text(size=16), axis.title=element_text(size=16), 
        legend.text=element_text(size=16),
        legend.key.size=unit(2,"lines"))& xlab(expression(beta)) & labs(col="", lwd="", lty="") 

ggsave("~/Dropbox/Generalized Count Splitting/JMLR-resubmit-sep-2023-v2/Figures/samp_split.png",
       width=12, height=4)


