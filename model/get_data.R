appr_url = "https://projects.fivethirtyeight.com/trump-approval-data/approval_polllist.csv"
elec_url = "https://projects.fivethirtyeight.com/polls-page/president_polls.csv"

get_approval = function() {
    pres_appr = suppressMessages(read_csv(appr_url)) %>%
        transmute(approval = approve/(approve + disapprove),
                  date = mdy(enddate)) %>%
        filter(date <= from_date) %>%
        arrange(date) %>%
        tail(20) %>%
        summarize(appr = mean(qlogis(approval))) %>%
        pull
}

get_elec_polls = function(write=F) {
    polls_raw = suppressWarnings(suppressMessages(read_csv(elec_url)))
    polls_d = polls_raw %>%
        group_by(question_id) %>%
        filter("Joseph R. Biden Jr." %in% candidate_name,
               "Donald Trump" %in% candidate_name,
               length(candidate_name) == 2,
               is.na(partisan), !internal) %>%
        ungroup %>%
        select(question_id, date1=start_date, date2=end_date, state, answer, pct,
                sample_size, population, pollster, fte_grade, tracking) %>%
        pivot_wider(names_from=answer, values_from=pct) %>%
        mutate(national = is.na(state),
               date = mdy(date1) + (mdy(date2) - mdy(date1))/2,
               day = n_days - ceiling(as.numeric(election_day - date) / 3) + 1,
               week = n_weeks - ceiling(as.numeric(election_day - date) / 21) + 1,
               dem = Biden / (Biden + Trump),
               gop = Trump / (Biden + Trump),
               sample = round(sample_size * (Biden + Trump)/100),
               logit_inflate = 1/dem + 1/(1 - dem),
               var_poll = logit_inflate^2 * dem * (1 - dem) / sample,
               firm = as.character(fct_lump(pollster, 50)),
               type_rv = population=="rv" | population == "v",
               type_lv = population=="lv",
               type_a = population=="a") %>%
        filter(date >= start_date, date <= from_date, !is.na(dem), firm != "Other") %>%
        select(state, national, date, day, week, firm, sample, 
               type_rv, type_lv, type_a, dem, var_poll)
    
    if (write) write_csv(polls_d, "data/polls/president_polls.csv")
    
    polls_d %>%
        left_join(state_abbr, by="state") %>%
        group_by(national) %>%
        mutate(abbr = if_else(national, "WA", abbr),
               date = as.character(date)) %>%
        select(abbr, day, date, week, everything(), -state) %>%
        mutate(state = match(abbr, state_abbr$abbr))    
}

get_gdp_est = function() {
    nowcast_url = str_c("https://www.newyorkfed.org/medialibrary/media/",
                        "research/policy/interactive/data/2020Q2.txt")
    nowcast = read_json(nowcast_url, simplifyVector=T) %>%
        as_tibble %>%
        filter(Year != "") %>%
        transmute(gdp_est = parse_number(Nowcast)) %>%
        pull
    
    list(gdp_est=tail(nowcast, 1)/100, gdp_sd=sd(nowcast/100))
}

get_state_prior_d = function() {
    state_abbr %>%
        mutate(year = 2020, inc = -1, two.term=0,
               pr_state = -(abbr=="FL") + (abbr=="DE"),
               vp_state = 0) %>%
        bind_rows(state_d, .) %>%
        group_by(state) %>%
        fill(id, region) %>%
        mutate(lpr_state = lag(pr_state, order_by=year),
               l2pr_state = lag(pr_state, 2, order_by=year),
               lvp_state = lag(vp_state, order_by=year),
               dem = qlogis(dem / (dem + gop)) - qlogis(natl_dem),
               ldem = lag(dem, order_by=year),
               l2dem = lag(dem, 2, order_by=year)) %>%
        filter(year=="2020") %>%
        select(year, state, region, ldem, l2dem, inc, two.term, pr_state, vp_state,
               lpr_state, lvp_state, l2pr_state) %>%
        left_join(state_abbr, by="state") %>%
        arrange(abbr)
}
