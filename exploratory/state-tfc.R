library(tidyverse)
library(rstanarm)

state_d = read_csv("data/historical/state_data_combined.csv")

d = state_d %>%
    group_by(state) %>%
    mutate(dem = qlogis(dem / (dem + gop)) - qlogis(natl_dem),
           lpr_state = lag(pr_state, order_by=year),
           l2pr_state = lag(pr_state, 2, order_by=year),
           lvp_state = lag(vp_state, order_by=year),
           ldem = lag(dem, order_by=year),
           l2dem = lag(dem, 2, order_by=year),
           inc = 2*(inc.party=="D") - 1,
           linc = lag(inc, order_by=year),
           ltt = lag(two.term, order_by=year),
           regn_yr = str_c(region, "_", year %% 100),
           q2.appr = qlogis(q2.appr)) %>%
    drop_na

cutoff_yr =  2016
d_test = filter(d, year >= cutoff_yr)
d_fit = filter(d, year < cutoff_yr)

lm1 = lm(dem ~ ldem + l2dem + pr_state + vp_state + lpr_state, data=d_fit)
lm2 = lm(dem ~ ldem + l2dem + pr_state + vp_state + lpr_state + l2pr_state + lvp_state, data=d_fit)
lm3 = lm(dem ~ ldem + l2dem + pr_state + vp_state + lpr_state + l2pr_state + lvp_state +
             inc + inc:two.term + (p_white + p_coll)*as.factor(year) + as.factor(region), data=d_fit)

m1 = stan_lmer(dem ~ ldem + l2dem + (1|region) + (1|year) +
                   pr_state + vp_state + lpr_state,
               data=d_fit, chains=1, warmup=500, iter=1700, prior=cauchy())
m2 = stan_lmer(dem ~ ldem + l2dem + (1|year:region) + inc + inc:two.term + 
                   pr_state + vp_state + lpr_state + l2pr_state + lvp_state,
               data=d_fit, chains=1, warmup=500, iter=1700, prior=cauchy())
m3 = stan_lmer(dem ~ ldem + l2dem + (1|region:year) + inc*ldem*two.term + 
                   pr_state + vp_state + lpr_state + l2pr_state + lvp_state +
                   (ldem||year),
               data=d_fit, chains=1, warmup=500, iter=1700, prior=cauchy())
m4 = stan_lmer(dem ~ ldem + l2dem + (ldem+l2dem|region) + (1|region:year) + 
                   inc*ldem*two.term +
                   pr_state + vp_state + lpr_state + l2pr_state + lvp_state,
               data=d_fit, chains=1, warmup=500, iter=1700, prior=cauchy())
m5 = stan_lmer(dem ~ ldem + l2dem + (ldem+l2dem|region) + (1|region:year) + 
                   inc*ldem*two.term + (region|year) +
                   pr_state + vp_state + lpr_state + l2pr_state + lvp_state,
               data=d_fit, chains=1, warmup=500, iter=1700, prior=cauchy())

print(m1, digits=3, detail=F)
print(m2, digits=3, detail=F)
print(m3, digits=3, detail=F)
print(m4, digits=3, detail=F)

compare_models(waic(m1), waic(m2), waic(m3), waic(m4))
compare_models(waic(m3), waic(m4), waic(m5))

ggplot(d_fit, aes(fitted(m5), resid(m5)/4, label=str_c(abbr, year %% 100))) +
    geom_text(size=3)

d %>%
    ungroup %>%
    mutate(fitted = colMeans(posterior_predict(m5, newdata=.))) %>%
    filter(year==2016) %>%
ggplot(aes(fitted, dem, label=str_c(abbr, year %% 100),
              color=sign(dem)==sign(fitted))) +
    guides(color=F) +
    geom_abline(slope=1) +
    xlim(NA, 1.2) +
    ylim(NA, 1.2) +
    geom_text(size=3)

spred_16 = posterior_predict(m5, newdata=d_test)
x = apply(spred_16, 2, function(x) x - mean(x))
xcov = (t(x) %*% x)/(nrow(x)-1)
xsd = sqrt(diag(xcov))
xcor = diag(1/xsd) %*% xcov %*% diag(1/xsd)
round(xcor[1:10, 1:10], 2)

heatmap(xcor, symm=T)
xcov = diag(xsd) %*% (0.8 + 0.2*diag(51)) %*% diag(xsd)

# FED GDP PREDICTIONS
# https://www.newyorkfed.org/medialibrary/media/research/policy/interactive/data/2020Q2.txt?nocache=1584541705932

# overall y vs fitted 
ggplot(d, aes(inc.share, est.share, label=lbl, color=correct)) + 
    geom_text(size=2) +
    geom_vline(xintercept=0.5) + geom_hline(yintercept=0.5)

# 2016 y vs fitted 
ggplot(d %>% filter(year==2016), 
       aes(inc.share, est.share, label=abbr, color=correct)) + 
    geom_text(size=2) +
    geom_vline(xintercept=0.5) + geom_hline(yintercept=0.5) +
    labs(title="2016 Election")

# EV: y vs fitted
ggplot(ev.winner, aes(pred.ev, inc.ev, label=year, color=(year >= cutoff.yr))) + 
    geom_text() +
    geom_vline(xintercept=270) + geom_hline(yintercept=270)

# incorrect y vs fitted 
ggplot(d %>% filter(!correct), 
       aes(est.share, inc.share, color=correct, label=lbl)) + 
    geom_hline(yintercept=0.5) +
    geom_text(size=1.5)



filter(state_d, year >= 2012) %>% 
    transmute(state = abbr, 
              year = year,
              dem_win = as.numeric(dem > gop), 
              margin = abs((dem - gop)/(dem + gop)))  %>% 
    pivot_wider(names_from=year, values_from=c(dem_win, margin)) %>%
    mutate_if(is.numeric, ~ round(., 4)) %>%
write_csv("docs/prev_results.csv", na="")








