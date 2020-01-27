library(tidyverse)
library(rstanarm)

d = read.csv("data/tfc.csv")

m1 = lm(vote ~ q2_growth + approval + incumbent, data=d)
m2 = stan_lm(vote ~ q2_growth + approval + incumbent, data=d, 
             prior=NULL, chains=1, iter=5000)

d.2020 = expand.grid(q2_growth=seq(-3, 3, 0.15), approval=seq(-20, 10, 0.75), 
                     incumbent = 0)

d.2020$pred = colMeans(posterior_predict(m2, d.2020) >= 50)

ggplot(d.2020) + 
    geom_tile(aes(x=approval, y=q2_growth, fill=pred)) +
    scale_fill_gradient2(low="#1040D0", mid="#EEEEEE", high="#CC1010", 
                         midpoint=0.5, position="right", limits=c(0.00, 1.00)) +
    labs(x="Net Approval, pct.", y="Q2 GDP Growth, pct.", fill="Pr(win)")

ggsave("heatmap.png", width=6, height=5)
