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

get_elec_polls = function(write=F, min_date=ymd("2020-03-01")) {
    keep_rule = function(d, key) {
        keep = (mdy_hm(d$created_at[1]) <= as_datetime(from_date)+86399)
        keep = keep & all(c("Joseph R. Biden Jr.", "Donald Trump") 
                     %in% d$candidate_name) & length(d$candidate_name) == 2
        keep = keep & (!d$internal[1] & 
            (!is.na(d$state) | str_detect(d$fte_grade[1], "[AB][+-]?"))) |
            (mdy(d$start_date[1]) <= ymd("2020-06-01") & d$population == "lv")
        
        d$keep = keep
        d
    }
    
    polls_raw = suppressWarnings(suppressMessages(read_csv(elec_url)))
    polls_d = polls_raw %>%
        filter(mdy(start_date) >= min_date) %>%
        group_by(question_id) %>%
        group_modify(keep_rule) %>%
        filter(keep) %>%
        ungroup %>%
        select(question_id, date1=start_date, date2=end_date, state, answer, pct,
                sample_size, population, pollster, fte_grade, tracking) %>%
        pivot_wider(names_from=answer, values_from=pct) %>%
        mutate(national = is.na(state),
               date = mdy(date1) + (mdy(date2) - mdy(date1))/2,
               day = n_days - ceiling(as.numeric(election_day - date) / 3),
               week = n_weeks - ceiling(as.numeric(election_day - date) / 21),
               dem = Biden / (Biden + Trump),
               gop = Trump / (Biden + Trump),
               sample_size = coalesce(sample_size, 400),
               sample = round(sample_size * (Biden + Trump)/100),
               logit_inflate = 1/dem + 1/(1 - dem),
               var_poll = logit_inflate^2 * dem * (1 - dem) / sample,
               firm = as.character(fct_lump(pollster, 50)),
               type_lv = population=="lv",
               type_a = population=="a",
               type_rv = !type_lv & !type_a) %>% #population=="rv" | population == "v") %>%
        filter(date >= start_date, date <= from_date, !is.na(dem), 
               !(firm == "Other" & date >= ymd("2020-08-01"))) %>%
        select(state, national, date, day, week, firm, sample, 
               type_rv, type_lv, type_a, dem, var_poll)
    
    if (write) write_csv(polls_d, "data/polls/president_polls.csv")
    
    polls_d %>%
        left_join(state_abbr, by="state") %>%
        left_join(state_regn, by="abbr") %>%
        group_by(national) %>%
        mutate(abbr = if_else(national, "WA", abbr),
               regn = if_else(national, "West", regn),
               date = as.character(date)) %>%
        select(abbr, regn, day, date, week, everything(), -state) %>%
        mutate(state = match(abbr, state_abbr$abbr)) %>%
        drop_na
}

get_gdp_est = function() {
    q2 = list(gdp_est=-0.314, gdp_sd=0.01)
    q3 = list(gdp_est=0.331, gdp_sd=0.02)

    #nowcast_url = str_c("https://www.newyorkfed.org/medialibrary/media/",
    #                    "research/policy/interactive/data/2020Q3.txt")
    #nowcast = read_json(nowcast_url, simplifyVector=T) %>%
    #    as_tibble %>%
    #    mutate(date = mdy(if_else(str_detect(Update, ", 2020"), 
    #                              Update, str_c(Update, ", 2020")))) %>%
    #    filter(date <= from_date, Year != "") %>%
    #    transmute(gdp_est = parse_number(Nowcast)) %>%
    #    pull
    
    #list(gdp_est=0.45*tail(nowcast, 1)/100 + 0.45*q2$gdp_est,
    #     gdp_sd=sqrt((0.2*0.15 + 0.3*sd(tail(nowcast, 5)/100))^2 + 0.5*q2$gdp_sd^2), 
    #     n=5)
    list(gdp_est=(1+q2$gdp_est)*(1+q3$gdp_est) - 1,
         gdp_sd=sqrt(q2$gdp_est^2 + q3$gdp_est^2),
         n=100)
}

get_state_prior_d = function() {
    state_abbr %>%
        mutate(year = 2020, inc = -1, two.term=0,
               pr_state = -(abbr=="NY") + (abbr=="DE"),
               vp_state = -(abbr=="IN") + (abbr=="CA")) %>%
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
