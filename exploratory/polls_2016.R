library(tidyverse)
library(lubridate)
library(tidybayes)
library(ggrepel)
library(usmap)

library(cmdstanr)
library(rstan)

# DATA LOADING  -------------------------------------------------------
election = ymd("2016-11-08")
start_date = ymd("2016-03-10")
n_days = ceiling(as.numeric(election - start_date) / 3)
n_weeks = ceiling(as.numeric(election - start_date) / 21)
to_date = function(day)  election - 3*(n_days - day + 1)

state_d = suppressMessages(read_csv("data/historical/state_data_combined.csv"))
states = distinct(select(state_d, state=abbr)) %>% pull 

run_date = ymd("2016-11-08")

state_reg = read_csv("data/historical/state_data_combined.csv") %>%
    select(state=abbr, regn=region) %>%
    distinct %>%
    left_join(distinct(select(state_d, name=state, state=abbr)), by="state") %>%
    select(state=name, regn)
state_reg$regn[state_reg$state=="Tennessee"] = "South"
state_reg$regn[state_reg$state=="West Virginia"] = "South"

polls_d = suppressMessages(read_csv("data/polls/president_general_polls_2016.csv")) %>%
    filter(type == "polls-only", 
           (state == "U.S." & str_sub(grade, 1, 1) == "A") 
           | (state != "U.S." & str_sub(grade, 1, 1) %in% c("A", "B", "C"))) %>%
    mutate(national = state == "U.S.",
           date = mdy(startdate) + (mdy(enddate) - mdy(startdate))/2,
           day = n_days - ceiling(as.numeric(election - date) / 3) + 1,
           week = n_weeks - ceiling(as.numeric(election - date) / 21) + 1,
           n_dem = round(samplesize * rawpoll_clinton/100),
           dem = rawpoll_clinton / (rawpoll_clinton + rawpoll_trump),
           gop = rawpoll_trump / (rawpoll_clinton + rawpoll_trump),
           sample = round(samplesize * (rawpoll_clinton + rawpoll_trump)/100),
           logit_inflate = 1/dem + 1/(1 - dem),
           var_poll = logit_inflate^2 * dem * (1 - dem) / sample,
           firm = as.character(fct_lump(pollster, 30)),
           type_rv = population=="rv" | population == "v",
           type_lv = population=="lv",
           type_a = population=="a") %>%
    filter(date >= start_date, !is.na(n_dem), firm != "Other") %>%
    filter(date <= run_date) %>%
    left_join(state_reg, by="state") %>%
    select(state, regn, national, date, day, week, firm, sample, 
           type_rv, type_lv, type_a, dem, n_dem, var_poll) 

# joined table + week stuff
wnum = (0:(n_days-1)) / 7
week_frac = if_else(floor(wnum)==wnum, 0, wnum - floor(wnum))
polls_d = polls_d %>%
    left_join(distinct(select(state_d, state, abbr)), by="state") %>%
    filter(!str_detect(state, "CD-")) %>%
    group_by(national) %>%
    mutate(abbr = if_else(national, "WA", abbr),
           regn = if_else(national, "West", regn),
           date = as.character(date)) %>%
    select(abbr, day, date, week, everything(), -state) %>%
    mutate(state = match(abbr, states))


# STAN DATA PREP  -------------------------------------------------------

poll_errors = read_rds("output/poll_errors.rdata")

model_d = compose_data(polls_d, .n_name = n_prefix("N"),
                       N_state = length(states),
                       D = n_days,
                       W = n_weeks,
                       D_W = ceiling(n_days/n_weeks),
                       week_frac,
                       week_day = floor(wnum) + 1,
                       prior_natl_mean = mean(pred_16), # from basic-tfc.R
                       prior_natl_sd = sd(pred_16),     # from   "
                       prior_state_mean = colMeans(spred_16), # from state-tfc.R
                       prior_state_cov = xcov,
                       #prior_state_cov = (0.08^2)*diag(N_state),
                       prior_rv_bias = 0.011,
                       prior_lv_bias = 0.0,
                       prior_a_bias = 0.02,
                       lv_rv_ratio = 8,
                       poll_errors,
                       )
#model_d$state = coalesce(model_d$state, 1) # fill natl polls (state=NA) with 1
model_d$prior_natl_poll_bias = 0
model_d$prior_all_state_poll_bias = 0
model_d$prior_regn_poll_bias = 0
str(model_d)


# MODEL FITTING -------------------------------------------------------

sm = cmdstan_model("stan/polls.stan", compiler_flags="CXXFLAGS += -Ofast")
#am = sm$optimize(data=model_d, iter=10000)
#write_stan_json(as.list(am$mle()), "stan/poll_mle.json")
sm$sample(data=model_d, num_chains=1, num_samples=1, num_warmup=500,
          #init="stan/poll_mle.json",
          adapt_delta=0.8, max_depth=15)$
    save_output_files("stan/", "polls_warmup", timestamp=F, random=F)
outf = read_file("stan/polls_warmup-1.csv")
stepsize = parse_number(str_extract(outf, "# Step size = [.0-9]+"))
imm = str_match(outf, "# Diagonal elements of inverse mass matrix:\n# ([ e\\-.,0-9]+)\n")[,2] %>% 
    str_split(", ", simplify=T) %>%
    parse_number

#sm = stan_model("stan/polls.stan")
#m = sampling(sm, data=model_d, chains=3, iter=1000, warmup=300, num_cores=4,
#             control=list(adapt_delta=0.95))

m = sm$sample(data=model_d, num_chains=3, num_samples=700, num_warmup=300,
              num_cores=4, adapt_delta=0.95, stepsize=stepsize, inv_metric=imm)
m = sm$sample(data=model_d, num_chains=3, num_samples=600, num_warmup=300,
              num_cores=4, adapt_delta=0.95)
m$save_output_files("stan/", "polls_warmup", timestamp=F, random=F)
m = sm$sample(data=model_d, num_chains=1, num_samples=600, num_warmup=300,
              num_cores=4, adapt_delta=0.95)

raw_draws = posterior::as_draws_df(m$draws())
uv_pars = c("sigma_natl", "sigma_state", "nonsamp_sd", "natl_error", 
            "all_state_error", "regn_error[1]", "regn_error[2]",
            "regn_error[3]", "regn_error[4]", "sd_firm", "bias_rv", "bias_lv")
raw_draws %>% select(uv_pars) %>%
    summarize_all(list(mean=mean, sd=sd)) %>% t

natl_draws = raw_draws %>%
    select(.chain, .iteration, .draw, starts_with("natl_dem")) %>%
    pivot_longer(cols=starts_with("natl_dem"), names_to="day", 
                 names_pattern="natl_dem\\[(.+)\\]", values_to="natl_dem") %>%
    mutate_if(is.character, as.numeric)

ggplot(natl_draws, aes(to_date(day), natl_dem)) + 
    stat_lineribbon() +
    geom_hline(yintercept=0.5111, lty="dashed") +
    geom_hline(yintercept=0.5, lty="solid") +
    scale_fill_brewer() +
    geom_line(aes(to_date(day), dem, group=firm, color=firm),
              data=filter(polls_d, national), 
              alpha=0.4, inherit.aes=F) + guides(color=F) +
    geom_text(aes(to_date(day), dem, label=abbreviate(firm, 3)),
              data=filter(polls_d, national), 
              alpha=0.5, size=2, inherit.aes=F)


state_draws = raw_draws %>%
    select(.chain, .iteration, .draw, starts_with("state_dem")) %>%
    pivot_longer(cols=starts_with("state_dem"), names_to=c("day", "state"),
                 names_pattern="state_dem\\[(.+),(.+)\\]", values_to="state_dem") %>%
    mutate(day = as.numeric(day),
           state_num = as.numeric(state),
           #state = state_abbr$abbr[as.numeric(state)])
           state = states[as.numeric(state)])

plot.state = "MI"
filter(state_draws, state==plot.state) %>%
ggplot(aes(to_date(day), state_dem)) + 
    stat_lineribbon() +
    geom_hline(yintercept=0.5, lty="solid") +
    scale_fill_brewer() +
    geom_line(aes(to_date(day), dem, group=firm, color=firm),
              data=filter(polls_d, abbr==plot.state, !national), 
              alpha=0.4, inherit.aes=F) + guides(color=F) +
    geom_text(aes(to_date(day), dem, label=abbreviate(firm, 3)), 
              data=filter(polls_d, abbr==plot.state, !national), 
              alpha=0.5, size=2, inherit.aes=F)

state_draws %>%
    group_by(state, day) %>%
    summarize(pr_win = mean(state_dem > 0.5)) %>%
    group_by(state) %>%
    filter(abs(pr_win[n_days] - 0.5) < 0.4) %>%
ggplot(aes(day, pr_win, color=state)) + 
    geom_line() +
    scale_color_viridis_d() +
    geom_text_repel(aes(label=state), data=~filter(.x, day==n_days),
                     segment.color="black", box.padding=0.1, nudge_x=6) +
    guides(color=F)
    

state_draws %>%
    filter(day == n_days) %>%
    group_by(state) %>%
    summarize(pr_win = mean(state_dem > 0.5)) %>%
ggplot(aes(pr_win, reorder(state, pr_win), label=state, color=pr_win)) +
    geom_text() + 
    guides(color=F) +
    scale_color_gradient2(low="red", high="blue", mid="#cccccc", midpoint=0.5)

state_draws %>%
    filter(day == n_days) %>%
    group_by(state) %>%
    summarize(pr_win = mean(state_dem > 0.5)) %>%
plot_usmap("states", value="pr_win", data=.) +
    scale_fill_gradient2(midpoint=0.5) +
    guides(fill=F)

state_ev = read_csv("data/historical/state_ev.csv") %>%
    left_join(distinct(select(state_d, state, abbr)), by="state") %>%
    select(state=abbr, ev=ev.2016)
evs = state_draws %>%
    filter(day == max(day)) %>%
    left_join(state_ev, by="state") %>%
    mutate(dem_ev = if_else(state_dem > 0.5, ev, 0)) %>%
    group_by(.draw) %>%
    summarize(dem_ev = sum(dem_ev))
ggplot(evs, aes(dem_ev, fill=dem_ev >= 270)) +
    geom_histogram(binwidth=1) +
    scale_fill_manual(values=c("red", "blue")) +
    guides(fill=F) + 
    geom_vline(xintercept=232, lty="dashed")
mean(evs$dem_ev >= 270)
mean(evs$dem_ev == 269)
mean(evs$dem_ev <= 232)

state_est_pct = state_draws %>%
    filter(day == n_days) %>%
    group_by(state) %>%
    summarize(dem_est = mean(state_dem)) 
filter(d, year==2016) %>%
    ungroup %>%
    transmute(state=abbr, dem_act = plogis(dem)) %>%
    left_join(state_est_pct, by="state") %>%
    mutate(err = dem_est - dem_act) %>%
    left_join(distinct(select(state_d, state=abbr, region)), by="state") %>%
ggplot(aes(dem_est, err, label=state, group=region, color=sign(dem_est-0.5))) +
    #geom_smooth(method=lm, se=F, color="black", size=0.5) +
    scale_color_gradient2() +
    geom_abline(slope=1) +
    geom_text(size=2) 

state_est_prob = state_draws %>%
    filter(day == n_days) %>%
    group_by(state) %>%
    summarize(pr_win = mean(state_dem > 0.5),
              dem_est  = mean(state_dem))
filter(d, year==2016) %>%
    ungroup %>%
    transmute(state=abbr, dem_act = plogis(dem)) %>%
    left_join(state_est_prob, by="state") %>%
    left_join(state_ev, by="state") %>%
    mutate(brier = ((dem_act > 0.5) - pr_win)^2,
           sqer = (dem_act - dem_est)^2) %>%
    summarize(brier_even = mean(brier),
              brier_ev = weighted.mean(brier, ev),
              rmse = sqrt(mean(sqer)))
 

state_draws %>%
    ungroup %>%
    filter(state == "VA" | state == "NC", day == n_days) %>%
    select(.draw, state, state_dem) %>%
    pivot_wider(names_from=state, values_from=state_dem) %>%
    rename(S1=2, S2=3) %>%
    summarize(mean(S1 > S2))
    #lm(PA ~ WI, data=.) %>%
    #summary()
ggplot(aes(TX, VA)) + geom_point(size=0.05, alpha=0.2) +
    coord_fixed() +
    geom_smooth(method=lm) +
    geom_abline(slope=1)

state_draws %>%
    filter(.draw == sample(1:1800, 1), day == n_days) %>%
    mutate(winner = round(state_dem)) %>%
plot_usmap("states", value="winner", data=.) +
    scale_fill_gradient2(midpoint=0.5) + #, limits=c(0.2, 0.8)) +
    guides(fill=F)

