library(tidyverse)
library(rstan)
library(magrittr)
library(loo)
library(bayesplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

d = read_csv("data/state_data_combined.csv")
d %<>% mutate(q2.pcinc = q2.inc / votes)

ev.winner = d %>% 
    group_by(year) %>%
    summarize(inc.ev = sum(ev[inc.share > 0.5]))
d %<>% left_join(ev.winner, by="year")

cutoff.yr = 2012
d.test = filter(d, year == cutoff.yr)
d.fit = filter(d, year < cutoff.yr)
T.fit = length(unique(d.fit$year))
T.test = length(unique(d.test$year))

years = unique(d$year)
states = unique(d$abbr)
to.int = function(x) as.integer(as.factor(x))

##############

sm = stan_model("state-result-model.stan")

model_fit = function(cutoff.yr, est=T, y_est=NULL) {
    d.fit = filter(d, year < cutoff.yr)
    T.fit = length(unique(d.fit$year))
    
    fit.data = with(d.fit, list(
        inc_party = 2*to.int(inc.party) - 3,
        year = to.int(year),
        state = to.int(state),
        X = model.matrix(~ 1 + q2.pcinc + q2.natl.gdp + q2.appr + two.term),
        y = if(is.null(y_est)) inc.share else y_est,
        
        estimate = est,
        prior_theta_scale = 5,
        prior_sigma_e_scale = 0.02,
        prior_sigma_t_scale = 0.1,
        prior_sigma_s_scale = 0.1,
        prior_sigma_s_location = 0.15,
        prior_eta_scale = 0.02
        ))
    fit.data$N = nrow(fit.data$X)
    fit.data$K = ncol(fit.data$X)
    fit.data$T = T.fit
    fit.data$S = 51
    
    fit = sampling(sm, data=fit.data, chains=1, iter=1000, warmup=600,
                   control=list(max_treedepth=15, adapt_delta=0.95))
    
    return(fit)
}

# FAKE DATA CHECK
est.fit = model_fit(2000, F)
fake.draws = rstan::extract(est.fit)
r.draw = 343
fake.fit = model_fit(2000, T, fake.draws$y_pred[r.draw,])

print(round(fake.draws$a_t[r.draw,], 3))
print(round(fake.draws$sigma_e[r.draw], 3))
print(round(fake.draws$sigma_t[r.draw], 3))
print(round(fake.draws$sigma_s[r.draw], 3))
print(round(fake.draws$beta[r.draw,], 3))
print(fake.fit, pars=c("theta", "log_lik", "y_pred", "a_s", "lin_pred"), include=F)
fake.a_s = colMeans(rstan::extract(fake.fit)$a_s)
qplot(fake.draws$a_s[r.draw,3,], fake.a_s[3,])

# MODEL CHECKING 
fit = model_fit(cutoff.yr)
draws = rstan::extract(fit)

print(fit, pars=c("theta", "log_lik", "y_pred", "a_s", "lin_pred", "A", "B"), include=F)

pairs(fit, pars=c("beta", "lp__"))
pairs(fit, pars=c("sigma_e", "sigma_t", "sigma_s", "lp__"))

# STATE EFFECTS PLOTTING

state.effects = colMeans(draws$a_s)
colnames(state.effects) = states
state.effects = as_tibble(state.effects) %>%
    mutate(year=years[1:nrow(.)], time.eft=colMeans(draws$a_t)) %>%
    gather(key="state", value="state.eft", -year, -time.eft)

ggplot(state.effects, aes(year, state.eft, group=state, label=state, 
                          alpha=(state %in% c("CA", "NE", "WV", "WA")))) + 
    geom_line() + guides(alpha=F) +
    geom_text(data=filter(state.effects, year==cutoff.yr-4), size=2, position="jitter")
ggplot(state.effects, aes(year, time.eft)) + geom_line()


#### predict new elections
model_predict = function(cutoff.yr, draws) {
    n = dim(draws$lp__)
    d.test = filter(d, year == cutoff.yr)
    T.test = 1
    T.fit = which.max(years == cutoff.yr) - 1
    test.data = with(d.test, list(
        inc_party = 2*to.int(inc.party) - 3,
        year = to.int(year),
        state = to.int(state),
        X = model.matrix(~ 1 + q2.pcinc + q2.natl.gdp + q2.appr + two.term),
        y = inc.share
        ))
    N = length(test.data$y)
    a_t = matrix(rnorm(n*T.test, 0, draws$sigma_t), n, T.test)
    a_s = array(rnorm(n*T.test*51, 0, 0.02), 
                dim=c(n, T.test, 51))
    a_s[,1,] = a_s[,1,] + draws$a_s[,T.fit,]
    #for (t in 2:T.test) {
    #    a_s[,t,] = a_s[,t,] + a_s[,t-1,]
    #}
    y_pred = array(dim=c(n, N))
    for (j in 1:N) {
        y_pred[,j] = with(test.data, rnorm(n, a_t[,year[j]] + draws$beta %*% X[j,]
                            + a_s[,year[j],state[j]]*inc_party[j], draws$sigma_e))
    }
    y_pred = boot::inv.logit(y_pred)
    
    return(y_pred)
}


y_pred = model_predict(cutoff.yr, draws)
evs = sweep(y_pred > 0.5, MARGIN=2, d.test$ev, `*`)
pop.votes = sweep(y_pred, MARGIN=2, d.test$votes, `*`)
distr = rowSums(evs)
pop.vote = rowSums(pop.votes)/sum(d.test$votes)
qplot(distr, bins=539, xlim=c(0,538)) 
qplot(distr, xlim=c(0,538)) 
qplot(pop.vote, bins=40)
median(distr)
mean(distr >= 270)
mean(pop.vote)

state.effects = rbind(colMeans(draws$a_s), colMeans(a_s))
colnames(state.effects) = states
state.effects = as_tibble(state.effects) %>%
    mutate(year=years[1:nrow(.)], 
           time.eft=c(colMeans(draws$a_t), colMeans(a_t))) %>%
    gather(key="state", value="state.eft", -year, -time.eft)

ggplot(state.effects, aes(year, state.eft, group=state, label=state, 
                          alpha=(state %in% c("CA", "NE", "WV", "WA")))) + 
    geom_line() + guides(alpha=F) +
    geom_text(data=filter(state.effects, year==cutoff.yr), size=2, position="jitter")
ggplot(state.effects, aes(year, time.eft)) + geom_line()

correct = sign(test.data$y - 0.5) == sign(colMeans(y_pred) - 0.5)
ggplot(NULL, aes(colMeans(y_pred), test.data$y, label=states, color=correct)) + 
    geom_text(size=3) +
    geom_vline(xintercept=0.5) + geom_hline(yintercept=0.5)

ggplot(NULL, aes(1:51, colMeans(y_pred > 0.5), label=states, color=correct)) +
    geom_text(size=3) +
    geom_hline(yintercept=0.5)

######################
#        TODO        #
######################
#  - prettify and automate model
#  - run back tests from 1984 or 1988 onward
#  - incorporate polling

# back tests
results = list()
for (t in 6:length(years)) {
    y = years[t]
    fit = model_fit(y)
    draws = rstan::extract(fit)
    y_pred = model_predict(y, draws)
    distr = rowSums(sweep(y_pred > 0.5, MARGIN=2, d.test$ev, `*`))
    results[[paste0("y_pred_", y)]] = y_pred
    results[[paste0("distr_", y)]] = distr
}

qplot(results$distr_1996)
qplot(results$distr_2000)
qplot(results$distr_2004)
qplot(results$distr_2008)
qplot(results$distr_2012)
qplot(results$distr_2016)

qplot(results$y_pred_1996[,10])
qplot(results$y_pred_2016[,10])

mean(results$distr_1996 >= 270)
mean(results$distr_2000 >= 270)
mean(results$distr_2004 >= 270)
mean(results$distr_2008 >= 270)
mean(results$distr_2012 >= 270)
mean(results$distr_2016 >= 270)
