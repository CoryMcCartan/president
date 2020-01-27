library(tidyverse)
library(lubridate)
library(magrittr)

state.returns = read_csv("data/state-returns/1976-2016-president.csv") %>%
    select(year, state, abbr=state_po, id=state_fips, party, 
           c.votes=candidatevotes, votes=totalvotes, writein) %>%
    mutate(party = recode(party, "democratic-farmer-labor"="democrat")) %>%
    filter(party == "democrat" | party == "republican", !writein) %>%
    spread(party, c.votes) %>%
    mutate(dem = democrat/votes,
           gop = republican/votes,
           margin = 0.5*(dem-gop)/(dem+gop)) %>%
    select(-writein, -democrat, -republican)

state.data = read_csv("data/state-data.csv") %>% select(-pop_est_2014)
state.returns %<>% left_join(state.data, by=c("state", "abbr"))

state.unemploy = read_csv("data/state_unemployment.csv") %>%
    transmute(id = str_sub(`Series ID`, 6, 7), 
           year = Year, 
           month = parse_integer(str_sub(Period, 2)),
           unemp = Value / 100) %>%
    filter(month >= 6, month <= 9, year %% 4 == 0) %>%
    spread(month, unemp) %>%
    rename(unemp.06 = `6`,
           unemp.07 = `7`,
           unemp.08 = `8`,
           unemp.09 = `9`)
state.returns %<>% left_join(state.unemploy, by=c("id", "year"))

state.gdp = read_csv("data/state_gdp.csv", na="(NA)", skip=4) %>%
    head(-7) %>%
    filter(LineCode == 1) %>%
    mutate(state = str_remove(GeoName, " \\*")) %>%
    select(-GeoName, -GeoFips, -Description, -LineCode) %>%
    gather(key="y.q", value="q2.inc", -state) %>%
    mutate(year = as.integer(str_sub(y.q, 1, 4)),
           quarter = as.integer(str_sub(y.q, 7))) %>%
    filter(year >= 1970, year %% 4 == 0, quarter == 2, 
        state != "United States") %>%
    select(state, year, q2.inc)
state.returns %<>% left_join(state.gdp, by=c("state", "year"))

pres.appr = read_csv("data/past_pres_approval.csv") %>%
    mutate(year = year(date),
           month = month(date),
           inc.party = if_else(party=="Democrat", "D", "R")) %>%
    filter(month == 6 | month == 5, year %% 4 == 0) %>%
    group_by(year, inc.party) %>%
    summarize(q2.appr = mean(approval)) %>%
    ungroup
state.returns %<>% left_join(pres.appr, by=c("year"))

tfc = read_csv("data/tfc.csv") %>%
    transmute(year=year, q2.natl.gdp = q2_growth/100, two.term = incumbent)
state.returns %<>% left_join(tfc, by=c("year"))

state.returns %<>% 
    mutate(inc.share = if_else(inc.party == "D", 0.5+margin, 0.5-margin))

state.ev = read_csv("data/state_ev.csv") %>%
    gather("year", "ev", -state) %>%
    mutate(year = parse_integer(str_sub(year, 4)))
state.returns %<>% left_join(state.ev, by=c("state", "year"))

write_csv(state.returns, "data/state_data_combined.csv")
