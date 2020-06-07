library(tidyverse)
library(lubridate)
library(rvest)
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

state.data = read_csv("data/historical/state-data.csv") %>% select(-pop_est_2014)
state.returns %<>% left_join(state.data, by=c("state", "abbr"))

state.unemploy = read_csv("data/historical/state_unemployment.csv") %>%
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

state.gdp = read_csv("data/historical/state_gdp.csv", na="(NA)", skip=4) %>%
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

state.growth = read_csv("data/historical/state_gdp.csv", na="(NA)", skip=4) %>%
    head(-7) %>%
    filter(LineCode == 1) %>%
    mutate(state = str_remove(GeoName, " \\*")) %>%
    select(-GeoName, -GeoFips, -Description, -LineCode) %>%
    gather(key="y.q", value="q2.inc", -state) %>%
    group_by(state) %>%
    mutate(year = as.integer(str_sub(y.q, 1, 4)),
           quarter = as.integer(str_sub(y.q, 7)),
           chg = q2.inc/lag(q2.inc, order_by=year) - 1,
           cycle = ceiling(year/4)*4) %>%
    filter(year > 1968, year %% 4 != 0) %>%
    group_by(state, cycle) %>%
    summarize(avg_growth = mean(chg)) %>%
    group_by(cycle) %>%
    mutate(natl_growth = avg_growth[which(state == "United States")]) %>%
    ungroup %>%
    filter(state != "United States") %>%
    transmute(state=state, year=cycle, avg_growth=avg_growth,
              diff_growth = avg_growth - natl_growth)
state.returns %<>% left_join(state.growth, by=c("state", "year"))


pres.appr = read_csv("data/polls/past_pres_approval.csv") %>%
    mutate(year = year(date),
           month = month(date),
           inc.party = if_else(party=="Democrat", "D", "R")) %>%
    filter(month == 6 | month == 5, year %% 4 == 0) %>%
    group_by(year, inc.party) %>%
    summarize(q2.appr = mean(approval)) %>%
    ungroup
state.returns %<>% left_join(pres.appr, by=c("year"))

tfc = read_csv("data/historical/tfc.csv") %>%
    transmute(year=year, q2.natl.gdp = q2_growth/100, two.term = incumbent)
state.returns %<>% left_join(tfc, by=c("year"))

state.returns %<>% 
    mutate(inc.share = if_else(inc.party == "D", 0.5+margin, 0.5-margin))

state.ev = read_csv("data/historical/state_ev.csv") %>%
    gather("year", "ev", -state) %>%
    mutate(year = parse_integer(str_sub(year, 4)))
state.returns %<>% left_join(state.ev, by=c("state", "year"))

state_d2 = haven::read_dta("data/historical/state_level_data.dta") %>%
    select(state=sstatename, abbr=sstate, p_evang, p_mormon, p_cath,
           p_metro=p_metroarea, p_black=pblack, p_hisp=phisp, p_coll=pbachelors) %>%
    mutate(state = if_else(state=="D.C.", "District of Columbia", state)) %>%
    mutate(p_white = 1-p_hisp-p_black)

w1972 = read_html("https://en.wikipedia.org/wiki/1972_United_States_presidential_election")
d1972 = w1972 %>% 
    html_nodes("table.wikitable") %>%
    `[[`(7) %>%
    html_table(header=F, fill=F) %>%
    slice(-1:-2, -54) %>%
    transmute(state = if_else(X1=="D.C.", "District of Columbia", X1),
           year = 1972,
           inc.party = "R",
           two.term = 0,
           dem = as.numeric(X7),
           dem = dem/(dem + as.numeric(X4)),
           gop = 1 - dem) %>%
    as_tibble %>%
    left_join(distinct(select(state.returns, state, abbr, id, region)), by="state")
w1968 = read_html("https://en.wikipedia.org/wiki/1968_United_States_presidential_election")
d1968 = w1968 %>% 
    html_nodes("table.wikitable") %>%
    `[[`(12) %>%
    html_table(header=F, fill=F) %>%
    slice(-1:-2, -54) %>%
    transmute(state = if_else(X1=="D.C.", "District of Columbia", X1),
           year = 1968,
           inc.party = "D",
           two.term = 1,
           dem = as.numeric(X7),
           dem = dem/(dem + as.numeric(X4)),
           gop = 1 - dem) %>%
    as_tibble %>%
    left_join(distinct(select(state.returns, state, abbr, id, region)), by="state")

home_states = read_csv("data/historical/home_states.csv")

ev_winner = state.returns %>% 
    group_by(year) %>%
    summarize(inc.ev = sum(ev[inc.share > 0.5]))

natl_d = read.csv("data/historical/tfc.csv") %>%
    select(-X,-X.1) %>%
    transmute(year=year,
              natl_inc = vote/100)

state.returns = state.returns %>%
    left_join(ev_winner, by="year") %>%
    left_join(state_d2, by=c("state", "abbr")) %>%
    bind_rows(d1968, d1972, .) %>%
    left_join(home_states, by="year") %>%
    left_join(natl_d, by="year") %>%
    group_by(state) %>%
    mutate(pr_state = (abbr == dem_pr_state) - (abbr == gop_pr_state),
           vp_state = (abbr == dem_vp_state) - (abbr == gop_vp_state),
           natl_dem = if_else(inc.party=="D", natl_inc, 1 - natl_inc),
           natl_gop = 1 - natl_dem,
           q2.pcinc = q2.inc/votes) %>%
    select(year, state, abbr, id, region, natl_dem, dem, gop, everything(), 
           -q2.inc, -dem_pr_state, -dem_vp_state, -gop_pr_state, -gop_vp_state,
           -natl_inc)
write_csv(state.returns, "data/historical/state_data_combined.csv")
