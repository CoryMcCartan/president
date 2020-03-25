library(tidyverse)
library(rstanarm)

d_tfc = read.csv("data/historical/tfc.csv") %>%
    select(-X,-X.1)

tfc_m = stan_lm(vote ~ q2_growth + approval + incumbent, data=d_tfc[-1,], 
             prior=R2(location=0.8), chains=1, iter=1000, 
             control=list(adapt_delta=0.999))

d.2020 = expand.grid(q2_growth=seq(-3, 3, 0.15), approval=seq(-20, 10, 0.75), 
                     incumbent = 0)

d.2020$pred = colMeans(posterior_predict(m, d.2020) >= 50)

ggplot(d.2020) + 
    geom_tile(aes(x=approval, y=q2_growth, fill=pred)) +
    scale_fill_gradient2(low="#1040D0", mid="#EEEEEE", high="#CC1010", 
                         midpoint=0.5, position="right", limits=c(0.00, 1.00)) +
    labs(x="Net Approval, pct.", y="Q2 GDP Growth, pct.", fill="Pr(win)")

ggsave("heatmap.png", width=6, height=5)

pred_16 = posterior_predict(m, d_tfc[1,])
mean(qlogis(pred_16/100))
sd(qlogis(pred_16/100))


