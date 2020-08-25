library(INLA)
library(tictoc)

polls_d = read_rds("output/polls_2020_08_22.rdata") %>%
    mutate(n_dem = round(dem*sample),
           type = case_when(type_rv ~ "rv", type_lv ~ "lv", T ~ "a"),
           id = 1:n(),
           abbr = if_else(national, "NATL", abbr),
           state = if_else(national, 52, state))

extra_d = tibble(
    day = (max(polls_d$day)+1):n_days,
    week = n_weeks - floor((n_days - day)/7),
    abbr = "NATL",
    state = 52,
    id = -1,
    firm = polls_d$firm[1],
    type = "lv",
    sample = 10000,
    n_dem = NA 
)

d = bind_rows(polls_d, extra_d)

{
tic()    
m40 = inla(n_dem ~ 1 + f(id, model="iid") + f(firm, model="iid") + 
               f(type, model="iid") + f(day, model="ar1") +
               f(abbr, group=week, model="iid", control.group=list(model="ar1")), 
           data=d, Ntrials=sample, family="binomial",
           control.family=list(control.link=list(model="logit")),
           control.compute=list(config = TRUE))
toc()
}


as_tibble(m40$summary.random$abbr) %>%
    mutate(day=rep(1:12, each=40),
           abbr=rep(ID[1:40], 12)) %>%
    filter(abbr=="MI") %>%
    pull(mean) %>% plot(type="l")

as_tibble(m40$summary.random$abbr) %>%
    rename(week=ID) %>%
    mutate(abbr=rep(ID[1:40], 12)) %>%
    filter(abbr=="MI") %>%
    pull(mean) %>% plot(type="l")