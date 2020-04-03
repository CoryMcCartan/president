#!/usr/bin/env Rscript

#####################################
#   U.S. PRESIDENTIAL MODEL 2020    #
#   CORY McCARTAN                   #
#   (c) 2020                        #
#####################################

library(optparse)

###
### Parse options
###

option_list = list(
    make_option("--dry", action="store_true", default=F,
                help="Dry run, results not saved."),
    make_option("--date", type="character", default=as.character(Sys.Date()),
                help="The date to estimate from."),
    make_option("--iter", type="integer", default=1800,
                help="Number of MCMC iterations for voter intent estimation,
                      not including warmup iterations."),
    make_option("--chains", type="integer", default=2,
                help="Number of MCMC chains for voter intent estimation."),
    make_option("--recompile", action="store_true", default=F,
                help="Force recompile of Stan models."),
    make_option("--model_dir", type="character", default="stan",
                help="The directory in which the models are stored"),
    make_option("--output_file", type="character", default="docs/estimate.json",
                help="The file to save estimates to."),
    make_option("--sims_file", type="character", default="docs/sims.json",
                help="The file to save state simulations to."),
    make_option("--history_file", type="character", default="docs/history.csv",
                help="The file to save model history to.")
)
opt = parse_args(OptionParser(option_list=option_list,
                              description="Forecast the 2020 U.S. presidential election."))
opt$iter = 3*ceiling(opt$iter/3)

library(tidyverse)
library(lubridate)
library(jsonlite)
library(tidybayes)
library(mvtnorm)
library(loo)
library(cmdstanr)
library(rstanarm)
options(mc.cores=4)

election_day = as.Date("2020-11-03")
start_date = ymd("2020-03-04")
from_date = as.Date(opt$date)
n_days = ceiling(as.numeric(election_day - start_date) / 3) + 1
n_weeks = ceiling(as.numeric(election_day - start_date) / 21) + 1
wnum = (0:(n_days-1)) / 7
week_frac = if_else(floor(wnum)==wnum, 0, wnum - floor(wnum))
to_date = function(day)  election_day - 3*(n_days - day)

polls_model_path = file.path(opt$model_dir, "polls")
natl_model_path = file.path(opt$model_dir, "natl-prior")
state_model_path = file.path(opt$model_dir, "state-prior")

state_abbr = read_rds("output/state_data_2016.rdata") %>%
    select(state, abbr) %>%
    filter(!str_detect(state, "CD-"))
state_regn = suppressMessages(read_csv("data/historical/state_data_combined.csv")) %>%
    select(abbr, regn=region) %>%
    distinct()
state_regn$regn[state_regn$abbr=="TN"] = "South"
state_regn$regn[state_regn$abbr=="WV"] = "South"

state_pop = suppressMessages(read_csv("data/state-returns/1976-2016-president.csv")) %>% 
    filter(year==2016) %>% 
    select(state=state_po, votes=totalvotes) %>%
    distinct()


###
### Download new polling data
###
cat("Downloading data.\n")
source("model/get_data.R")

# presidential approval, general election polls, and distribution of past polling errors
pres_appr = get_approval() # mean pct. pt. gap appr-disappr
polls_d = get_elec_polls()

poll_errors = read_rds("output/poll_errors.rdata")
poll_errors$prior_natl_poll_bias = 0
poll_errors$prior_all_state_poll_bias = 0
poll_errors$prior_regn_poll_bias = 0

# Q2 GDP forecast
gdp_growth = get_gdp_est()
gdp_samp = rt(500, gdp_growth$n-1)*gdp_growth$gdp_sd + gdp_growth$gdp_est

# past state results and other covariates
state_d = suppressMessages(read_csv("data/historical/state_data_combined.csv"))
natl_prior_d2020 = tibble(year=2020, q2_gdp=gdp_samp, q2_appr=pres_appr, two_term=0)
state_prior_d2020 = get_state_prior_d()

###
### Fitting prior model
###
cat("Building prior.\n")
source("model/get_models.R")

natl_model = get_natl_prior_m(natl_model_path, opt$recompile)
natl_prior_pred = -as.numeric(posterior_predict(natl_model, newdata=natl_prior_d2020))

state_model = get_state_prior_m(state_model_path, state_d, opt$recompile)
state_prior_pred = posterior_predict(state_model, newdata=state_prior_d2020)
state_prior_mean = colMeans(state_prior_pred)
state_x = apply(state_prior_pred, 2, function(x) x - mean(x))
state_prior_cov = (t(state_x) %*% state_x) / (nrow(state_x) - 1)
state_sd = sqrt(diag(state_prior_cov))
state_prior_corr = diag(1/state_sd) %*% state_prior_cov %*% diag(1/state_sd)

###
### Fitting main model
###
cat("Fitting main model.\n")

model_d = compose_data(polls_d, .n_name = n_prefix("N"),
                       N_state = nrow(state_abbr),
                       D = n_days,
                       W = n_weeks,
                       D_W = ceiling(n_days/n_weeks),
                       week_frac,
                       week_day = floor(wnum) + 1,
                       prior_natl_mean = mean(natl_prior_pred), 
                       prior_natl_sd = sqrt(2)*sd(natl_prior_pred),
                       prior_state_mean = state_prior_mean, 
                       prior_state_cov = state_prior_cov,
                       prior_rv_bias = 0.011,
                       prior_lv_bias = 0.0,
                       prior_a_bias = 0.02,
                       lv_rv_ratio = 8,
                       poll_errors,
)

polls_model = get_polls_m(polls_model_path, opt$recompile)

# TODO incorporate inv_metric stuff
fit_polls = polls_model$sample(data=model_d, num_chains=3, num_samples=opt$iter/3, 
                               num_warmup=300, num_cores=4, adapt_delta=0.97, 
                               stepsize=0.015)

draws = posterior::as_draws_df(fit_polls$draws())

###
### Predictions
###

natl_draws = draws %>%
    select(.chain, .iteration, .draw, starts_with("natl_dem")) %>%
    pivot_longer(cols=starts_with("natl_dem"), names_to="day", 
                 names_pattern="natl_dem\\[(.+)\\]", values_to="natl_dem") %>%
    mutate_if(is.character, as.numeric)

natl_final = filter(natl_draws, day==n_days) %>% pull(natl_dem)

state_draws = draws %>%
    select(.chain, .iteration, .draw, starts_with("state_dem")) %>%
    pivot_longer(cols=starts_with("state_dem"), names_to=c("day", "state"),
                 names_pattern="state_dem\\[(.+),(.+)\\]", values_to="state_dem") %>%
    mutate(day = as.numeric(day),
           state_num = as.numeric(state),
           state = state_abbr$abbr[state_num])

state_ev = suppressMessages(read_csv("data/historical/state_ev.csv")) %>%
    left_join(state_abbr, by="state") %>%
    select(state=abbr, ev=ev.2016)
evs = state_draws %>%
    filter(day == max(day)) %>%
    left_join(state_ev, by="state") %>%
    mutate(dem_ev = if_else(state_dem > 0.5, ev, 0)) %>%
    group_by(.draw) %>%
    summarize(dem_ev = sum(dem_ev)) %>%
    pull

# IMPORTANCE SAMPLING

{
#tgt_param = draws$`mu_natl[83]`
#lprior_now = with(model_d, dnorm(tgt_param, prior_natl_mean, prior_natl_sd, log=T))
#lprior_tgt = with(model_d, dt((tgt_param - prior_natl_mean)/prior_natl_sd, 3, log=T))
    
tgt_param = draws$bias_rv
lprior_now = with(draws, dnorm(bias_rv, model_d$prior_rv_bias, 0.04, log=T) + 
                      dnorm(bias_lv, 0, 0.04/model_d$lv_rv_ratio, log=T))
lprior_tgt = with(draws, dnorm(bias_rv, -0.011, 0.01, log=T) + 
                      dnorm(bias_lv, 0, 5e-3, log=T))
    
lratio = lprior_tgt - lprior_now
qplot(lratio)
r_eff = relative_eff(exp(-lratio), draws$.chain)

psis_wgt = psis(lratio, r_eff=r_eff) 
print(pareto_k_values(psis_wgt))
is_wgt = weights(psis_wgt, log=F)

print(mean(natl_final))
print(weighted.mean(natl_final, is_wgt))

print(mean(evs >= 270))
print(weighted.mean(evs >= 270, is_wgt))

print(model_d$prior_rv_bias)
print(mean(draws$bias_rv))
print(weighted.mean(draws$bias_rv, is_wgt))

qplot(natl_final, is_wgt) 
qplot(evs, is_wgt) 
}


{
tgt_param =  draws %>%
    select(starts_with("mu_state[13,")) %>%
    as.matrix

rho = median(state_prior_corr[!diag(51)])
#cov_test = diag(state_sd) %*% (rho + diag(51)*(1-rho)) %*% diag(state_sd)
corr_test = rlang::duplicate(state_prior_corr)
corr_test[!diag(51)] = pmin(corr_test[!diag(51)] * 1.5, 0.99)
cov_test = diag(state_sd) %*% corr_test %*% diag(state_sd)
#cov_test = rlang::duplicate(state_prior_cov)
#cov_test[!diag(51)] = cov_test[!diag(51)] * 1.3
    
lprior_now = with(model_d, dmvnorm(tgt_param, prior_state_mean, prior_state_cov, log=T))
lprior_tgt = with(model_d, dmvnorm(tgt_param, prior_state_mean, cov_test, log=T))
lratio = lprior_tgt - lprior_now
qplot(lratio)
r_eff = relative_eff(exp(-lratio), draws$.chain)

psis_wgt = psis(lratio, r_eff=r_eff) 
print(pareto_k_values(psis_wgt))
is_wgt = weights(psis_wgt, log=F)

print(mean(natl_final))
print(weighted.mean(natl_final, is_wgt))
print(weighted.mean(natl_final, exp(lratio)))

print(mean(evs >= 270))
print(weighted.mean(evs >= 270, is_wgt))

qplot(evs, is_wgt)
hist = map_dbl(0:538, ~ mean(evs == .))
hist_wgt = map_dbl(0:538, ~ sum(is_wgt[evs == .]))
qplot(0:538, hist, fill=(0:538 >= 270), geom="col") + guides(fill=F)
qplot(0:538, hist_wgt, fill=(0:538 >= 270), geom="col") + guides(fill=F)
}





