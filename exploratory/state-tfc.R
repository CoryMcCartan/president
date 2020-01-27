library(tidyverse)
library(rstanarm)
library(lme4)
library(magrittr)

d = read_csv("data/state_data_combined.csv")
d %<>% mutate(q2.pcinc = q2.inc / votes)

ev.winner = d %>% 
    group_by(year) %>%
    summarize(inc.ev = sum(ev[inc.share > 0.5]))
d %<>% left_join(ev.winner, by="year")

cutoff.yr =  2008
d.test = filter(d, year >= cutoff.yr)
d.fit = filter(d, year <cutoff.yr)

m1 = lmer(inc.share ~ (0 + inc.party | abbr) 
          + (0 + inc.party | region) + (1 | year)
          + q2.pcinc + q2.natl.gdp + q2.appr + two.term, data=d.fit)

m3 = lmer(inc.share ~ (0 + inc.party | abbr) 
          + (1 | year) + q2.pcinc + q2.natl.gdp + q2.appr + two.term, data=d.fit)


m1b = stan_lmer(inc.share ~ (0 + inc.party | abbr) 
          + (0 + inc.party | region) + (1 | year)
          + q2.pcinc + q2.natl.gdp + q2.appr + two.term, data=d.fit,
          chains=1, iter=1000)



yh.b = posterior_predict(m1b, d)
evs = sweep(yh.b > 0.5, MARGIN=2, d$ev, `*`)
distr.2016 = rowSums(evs[,511:561])
distr.2012 = rowSums(evs[,460:510])
distr.2008 = rowSums(evs[,409:459])


d %<>% mutate(est.share = predict(m3, d, allow.new.levels=T),
              correct = sign(est.share - 0.5) == sign(inc.share - 0.5),
              lbl = paste(abbr, year))

ev.winner = d %>% 
    group_by(year) %>%
    summarize(pred.ev = sum(ev[est.share > 0.5])) %>%
    right_join(ev.winner, by="year")

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
