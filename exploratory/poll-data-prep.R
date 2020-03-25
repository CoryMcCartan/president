library(tidyverse)

### BAYES MACHINERY

# PRIOR SETUP
mu_0 = 0
n_mu = 4
prec_0 = 100
n_prec = 2

# plot priors
plot(function(x) dnorm(x, mu_0, sqrt(1/n_mu/prec_0)), xlim=c(-0.5, 0.5))
qplot(rgamma(1e4, n_prec/2, rate=n_prec/(2*prec_0))^(-0.5), bins=80) + xlim(0, 1)

plot(function(x) dt((x - -0.00299)/0.0788, 6), xlim=c(-0.5, 0.5))

round((7.8*rt(10, 6) - 0.299)/4, 1)
round((2.13*rt(10, 231) - 1.36)/4, 1)


post_error = function(xbar, s, n, ...) {
    mu = (n_mu*mu_0 + n*xbar)/(n_mu + n)
    n_mu2 = n_mu + n 
    alpha = (n_prec + n)/2
    beta = (n_prec/prec_0 + n*s^2 + n_mu*n*(xbar - mu_0)^2/(n_mu+n))/2
    list(nu = n_mu + n, mu = mu, sigma = sqrt(beta/alpha), sd=s)
}

#### POLL DATA ########

raw_polls = read_csv("data/polls/raw-polls.csv") %>%
    filter(type_simple == "Pres-G") %>%
    mutate(national = location == "US",
           wgt = 100^2 * samplesize / (cand1_pct*(100-cand1_pct)),
           dem_poll = cand1_pct / (cand1_pct + cand2_pct),
           dem_act = cand1_actual / (cand1_actual + cand2_actual)) %>%
    select(year, state=location, national, samplesize, wgt, dem_poll, dem_act)

polls = raw_polls %>%
    group_by(year, state, national) %>%
    summarize(dem_poll = weighted.mean(dem_poll, wgt),
           dem_act = mean(dem_act))


## NATL ERROR

# avg. natl error and bias
natl_error = polls %>%
    filter(national) %>%
    ungroup %>%
    mutate(logit_error = qlogis(dem_poll) - qlogis(dem_act)) %>%
    summarize(xbar = mean(logit_error), s = sd(logit_error), n = n()) %>%
    pmap_dfr(post_error)

# elec-by-ele national error
natl_error_elec = polls %>%
    filter(national) %>%
    group_by(year) %>%
    transmute(natl_error = qlogis(dem_poll) - qlogis(dem_act))

state_error_elec = polls %>%
    ungroup %>%
    filter(!national) %>%
    left_join(natl_error_elec, by="year") %>%
    mutate(logit_error = qlogis(dem_poll) - qlogis(dem_act) - natl_error) %>%
    group_by(year) %>%
    summarize(state_error = mean(logit_error))

all_state_errors = state_error_elec %>%
    summarize(xbar = mean(state_error), s = sd(state_error), n = n()) %>%
    pmap_dfr(post_error)
    
states_errors = polls %>%
    ungroup %>%
    filter(!national) %>%
    left_join(state_error_elec, by="year") %>%
    left_join(natl_error_elec, by="year") %>%
    mutate(logit_error = qlogis(dem_poll) - qlogis(dem_act) - state_error - natl_error) %>%
    summarize(xbar = mean(logit_error), s = sd(logit_error), n = n()) %>%
    pmap_dfr(post_error)

infl_factor = sd(rt(1e6, natl_error$nu))
poll_errors = list(prior_natl_poll_bias = natl_error$mu,
                   prior_natl_poll_error = infl_factor * natl_error$sigma,
                   prior_all_state_poll_bias = all_state_errors$mu,
                   prior_all_state_poll_error = infl_factor * all_state_errors$sigma,
                   prior_states_poll_error = states_errors$sigma)
write_rds(poll_errors, "output/poll_errors.rdata")

cm = polls %>%
    filter(!national) %>%
    left_join(natl_error_elec, by="year") %>%
    mutate(logit_error = qlogis(dem_poll) - qlogis(dem_act)) %>%
    select(year, state, logit_error) %>%
    pivot_wider(names_from=state, values_from=logit_error) %>%
    ungroup %>%
    select(-year) %>%
    as.matrix %>%
    cor(use="pairwise.complete.obs") 
heatmap(cm[-51,-51], symm=T)

polls %>%
    ungroup %>%
    filter(!national) %>%
    mutate(logit_error = qlogis(dem_poll) - qlogis(dem_act)) %>%
    group_by(year) %>%
    summarize(xbar = mean(logit_error), s = sd(logit_error), n=n())
polls %>%
    ungroup %>%
    filter(!national) %>%
    mutate(logit_error = qlogis(dem_poll) - qlogis(dem_act)) %>%
    summarize(xbar = mean(logit_error), s = sd(logit_error), n=n())

polls %>%
    ungroup %>%
    filter(!national) %>%
    mutate(logit_error = qlogis(dem_poll) - qlogis(dem_act)) %>%
aov(logit_error ~ as.factor(year) + as.factor(state), data=.) %>%
    anova

state.d = read_rds("output/state_data_2016.rdata") %>%
    select(state=abbr, everything(), -state, -pop) %>%
    mutate(votes = log(votes))

library(usmap)
polls %>% filter(year==2016) %>% mutate(error = dem_poll-dem_act) %>%
    plot_usmap("states", values="error", data=.) +
    scale_fill_gradient2() +
    guides(fill=F)
