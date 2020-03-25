library(tidyverse)
library(censusapi)

state.d = read_csv("data/historical/state_data_combined.csv") %>%
    filter(year == 2012) %>%
    bind_rows(
        tibble(state="Nebraska CD-2", abbr="NE-2", dem=0.4570, gop=0.5285),
        tibble(state="Maine CD-2", abbr="ME-2", dem=0.5294, gop=0.4438)
    ) %>%
    mutate(dem_shift = qlogis(dem/(dem+gop)) - qlogis(0.5198)) %>%
    select(state, abbr, votes, dem_shift) %>%
    bind_rows(tibble(state="U.S.", abbr="US", dem_shift=0)) %>%
    arrange(abbr)

census.vars = c("NAME", "B01001_001E", "B02001_002E", 
                "B16010_041E", "B19013_001E")
census_year = 2014
state.demg = getCensus("acs/acs1", census_year, region="state:*", vars=census.vars)
nebr2.demg = getCensus("acs/acs1", census_year, region="congressional district:02", 
                       regionin="state:31", vars=census.vars) %>%
    mutate(NAME="Nebraska CD-2") %>%
    select(-congressional_district)
maine2.demg = getCensus("acs/acs1", census_year, region="congressional district:02", 
                        regionin="state:23", vars=census.vars) %>%
    mutate(NAME="Maine CD-2") %>%
    select(-congressional_district)
state.demg = bind_rows(state.demg, nebr2.demg, maine2.demg) %>%
    transmute(state = NAME,
              pop = B01001_001E,
              pct_white = scale(B02001_002E / pop, scale=F),
              pct_coll = scale(B16010_041E / pop, scale=F),
              medinc_k = scale(B19013_001E))

state.d = inner_join(state.d, state.demg, by="state")

write_rds(state.d, "output/state_data_2016.rdata", compress="gz")



# SYNTHETIC ELECTORATES

raw = read_csv("~/Documents/Analyses/census/acs_p.csv")
states = read_csv("~/Documents/Analyses/census/states.csv")

acs = raw %>%
    mutate(AGEP = as.integer(AGEP)) %>%
    filter(AGEP >= 18) %>%
    filter(CIT != 5) %>%
    left_join(states, by=c("ST"="id")) %>%
    transmute(
        abbr = abbr,
        age = case_when(
            AGEP < 35 ~ "18-34",
            AGEP < 49 ~ "35-49",
            AGEP < 65 ~ "50-64",
            T ~ "65+"
        ),
        educ = case_when(
            SCHL <= 16 ~ "hs",
            SCHL <= 20 ~ "somecoll",
            T ~ "coll"
        ),
        inc = as.numeric(PINCP),
        inc = case_when(
            inc <= 30e3 ~ "$",
            inc <= 60e3 ~ "$$",
            inc <= 100e3 ~ "$$$",
            inc <= 200e3 ~ "$$$$",
            T ~ "$$$$$"
        ),
        hisp = HISP != "01",
        race = RAC1P %>% as.factor %>% fct_collapse(white = "1", black = "2",
            other = c("3", "4", "5", "6", "7", "8", "9")) %>%
            as.character %>%
            if_else(hisp, "hisp", .),
        weight = as.integer(PWGTP)
    ) %>% 
    select(-hisp)

demg_matrix = acs %>% 
    drop_na %>%
    group_by(abbr, age, educ, inc, race) %>%
    summarize(weight = sum(weight)) %>%
    group_by(abbr) %>%
    mutate(weight = weight / sum(weight)) %>%
    pivot_wider(names_from=abbr, values_from=weight, values_fill=list(weight=0)) %>%
    select(-age, -educ, -inc, -race) %>%
    as.matrix %>% t

state_cor = cov2cor(demg_matrix %*% t(demg_matrix))

x = c("AL", "CA", "FL", "MN", "NC", "NM", "RI", "WI")
round(state_cor[x, x], 2)

write_rds(state_cor, "output/state_demg_corr.rdata", compress="gz")
