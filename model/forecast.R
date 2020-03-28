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
    make_option("--iter", type="integer", default=2000,
                help="Number of MCMC iterations for voter intent estimation,
                      including 500 warmup iterations."),
    make_option("--chains", type="integer", default=2,
                help="Number of MCMC chains for voter intent estimation."),
    make_option("--recompile", action="store_true", default=F,
                help="Force recompile of Stan models."),
    make_option("--model_dir", type="character", default="stan",
                help="The directory in which the models are stored"),
    make_option("--output_file", type="character", default="docs/estimate.json",
                help="The file to save estimates to."),
    make_option("--history_file", type="character", default="docs/history.csv",
                help="The file to save model history to.")
)
opt = parse_args(OptionParser(option_list=option_list,
                              description="Forecast the 2020 U.S. presidential election."))

suppressMessages(library(tibble))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(stringr))
suppressMessages(library(forcats))
suppressMessages(library(readr))
suppressMessages(library(lubridate))
suppressMessages(library(glue))
suppressMessages(library(tidybayes))
suppressMessages(library(jsonlite))

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


###
### Download new polling data
###
cat("Downloading data.\n")
source("model/get_data.R")

# presidential approval, general election polls, and distribution of past polling errors
pres_appr = get_approval() # mean pct. pt. gap appr-disappr
old_polls = suppressMessages(read_csv("docs/polls.csv", col_types="cilcd"))
polls_d = get_elec_polls()
polls_d %>%
    select(date, state, national, firm, dem) %>%
    write_csv("docs/polls.csv")
if (all.equal(old_polls, select(polls_d, date, state, national, firm, dem))) {
    cat("No new polls.\n")
    system("osascript -e beep"); system("osascript -e beep")
    Sys.sleep(10)
}
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
suppressMessages(library(cmdstanr))
suppressMessages(library(rstanarm))
options(mc.cores=4)
source("model/get_models.R")

natl_model = get_natl_prior_m(natl_model_path, opt$recompile)
natl_prior_pred = -as.numeric(posterior_predict(natl_model, newdata=natl_prior_d2020))

state_model = get_state_prior_m(state_model_path, state_d, opt$recompile)
state_prior_pred = posterior_predict(state_model, newdata=state_prior_d2020)
state_prior_mean = colMeans(state_prior_pred)
state_x = apply(state_prior_pred, 2, function(x) x - mean(x))
state_prior_cov = (t(state_x) %*% state_x) / (nrow(state_x) - 1)
#state_sd = sqrt(diag(state_prior_cov))
#state_prior_cov = diag(state_sd) %*% (0.9 + 0.1*diag(51)) %*% diag(state_sd)

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
fit_polls = polls_model$sample(data=model_d, num_chains=3, num_samples=700, 
                               num_warmup=300, num_cores=4, adapt_delta=0.97, 
                               stepsize=0.015)

raw_draws = posterior::as_draws_df(fit_polls$draws())

###
### Output predictions
###

natl_draws = raw_draws %>%
    select(.chain, .iteration, .draw, starts_with("natl_dem")) %>%
    pivot_longer(cols=starts_with("natl_dem"), names_to="day", 
                 names_pattern="natl_dem\\[(.+)\\]", values_to="natl_dem") %>%
    mutate_if(is.character, as.numeric)

natl_final = filter(natl_draws, day==n_days) %>% pull(natl_dem)

state_draws = raw_draws %>%
    select(.chain, .iteration, .draw, starts_with("state_dem")) %>%
    pivot_longer(cols=starts_with("state_dem"), names_to=c("day", "state"),
                 names_pattern="state_dem\\[(.+),(.+)\\]", values_to="state_dem") %>%
    mutate(day = as.numeric(day),
           state_num = as.numeric(state),
           state = state_abbr$abbr[state_num])

firm_ct = polls_d %>%
    group_by(firm) %>%
    summarize(n=n())
firm_ids = tibble(firm=polls_d$firm, id=model_d$firm) %>%
    distinct()
firms = raw_draws %>%
    select(.chain, .iteration, .draw, starts_with("house_effects")) %>%
    pivot_longer(cols=starts_with("house_effects"), names_to="id", 
             names_pattern="house_effects\\[(.+)\\]", values_to="effect") %>%
    mutate_if(is.character, as.numeric) %>%
    group_by(id) %>%
    summarize(effect = median(effect)) %>%
    left_join(firm_ids, by="id") %>%
    left_join(firm_ct, by="firm") %>%
    select(firm, n, effect)

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


entry = tibble(
    date = from_date,
    ev_exp = median(evs),
    ev_min = min(evs),
    ev_q05 = quantile(evs, 0.05),
    ev_q25 = quantile(evs, 0.25),
    ev_q75 = quantile(evs, 0.75),
    ev_q95 = quantile(evs, 0.95),
    ev_max = max(evs),
    natl_exp = median(natl_final),
    natl_q05 = quantile(natl_final, 0.05),
    natl_q25 = quantile(natl_final, 0.25),
    natl_q75 = quantile(natl_final, 0.75),
    natl_q95 = quantile(natl_final, 0.95),
    prob = mean(evs >= 270),
    prob_pop = mean(natl_final > 0.5),
    prob_tie = mean(evs == 269),
    prob_pop_ev_DR = mean(natl_final > 0.5 & evs <= 268),
    prob_pop_ev_RD = mean(natl_final < 0.5 & evs >= 270)
)
entry = mutate_if(entry, is.numeric, ~ round(., 4))

natl_intent_tbl = natl_draws %>%
    group_by(day) %>%
    summarize(natl_exp = median(natl_dem),
              natl_q05 = quantile(natl_dem, 0.05),
              natl_q25 = quantile(natl_dem, 0.25),
              natl_q75 = quantile(natl_dem, 0.75),
              natl_q95 = quantile(natl_dem, 0.95)) %>%
    mutate(date = to_date(day))

state_probs = state_draws %>%
    filter(day == n_days) %>%
    group_by(state) %>%
    summarize(prob = mean(state_dem > 0.5)) %>%
    left_join(state_abbr, by=c("state"="abbr")) %>%
    left_join(state_ev, by="state") %>%
    rename(state_name = state.y)

output = append(as.list(entry), list(
    time = Sys.time(),
    n_polls = nrow(polls_d),
    natl_intent = natl_intent_tbl,
    firm_effects = firms,
    states = state_probs,
    hist = map_int(0:538, ~ sum(evs == .))
))

with(output, cat(glue("
    ===========================================
     2020 U.S. Presidential Forecast
     {as.character(Sys.Date(), format='%B %d, %Y')}
    -------------------------------------------
     Forecast from: {as.character(from_date, format='%B %d, %Y')}
     {round((election_day - from_date)/7)} weeks until the election.
     {n_polls} polls.

     Dem. share of two-party vote:   {round(100*natl_exp, 1)}%
     Median EV estimate:             {round(ev_exp)}
     Estimated EV range:             {round(ev_q05)} - {round(ev_q95)}
     Probability of winning:         {round(100*prob)}%
    ===========================================
    ")))
cat("\n\n")

if (opt$dry) quit("no")

if (file.exists(opt$history_file)) {
    history = bind_rows(
        suppressMessages(read_csv(opt$history_file)),
        entry
    ) %>%
        group_by(date) %>%
        slice(n()) %>%
        ungroup %>%
        arrange(date)
} else {
    history = entry
}
write_csv(history, opt$history_file, na="")

# only save full output if current run
if (from_date != Sys.Date()) quit("no")
write_json(output, opt$output_file, auto_unbox=T, digits=7)
