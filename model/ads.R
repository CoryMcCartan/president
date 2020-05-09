#!/user/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(jsonlite))
suppressMessages(library(usmap))
suppressMessages(library(httr))

# setup basic data
state_abbr = read_rds("output/state_data_2016.rdata") %>%
    select(state, abbr) %>%
    filter(!str_detect(state, "CD-"))
state_abbr$state[state_abbr$abbr=="DC"] = "Washington, District of Columbia"

state_pop = suppressMessages(read_csv("data/state-returns/1976-2016-president.csv")) %>% 
    filter(year==2016) %>% 
    select(state=state_po, votes=totalvotes) %>%
    distinct()

states_2016 = read_csv("data/historical/state_data_combined.csv") %>% 
    filter(year==2016) %>% 
    select(abbr, dem_16=dem)

page_ids = c(biden="7860876103", trump="153080620724")

get_token = function() {
    app_id = "3080333892192928"
    app_secret = Sys.getenv("FB_PRES_2020_SECRET")
    fb_app = oauth_app(appname="facebook", key=app_id, secret=app_secret)
    resp = oauth2.0_token(oauth_endpoints("facebook"), fb_app, scope="public_profile",
                   type="application/x-www-form-urlencoded", cache=TRUE)
    fromJSON(names(resp$credentials))$access_token
}

get_url = function(cand_id, num=1000)  {
    str_glue("https://graph.facebook.com/v6.0/ads_archive?",
                   "access_token={get_token()}&",
                   "ad_type=POLITICAL_AND_ISSUE_ADS&",
                   "search_page_ids=[{paste0(cand_id, collapse=',')}]&",
                   "ad_reached_countries=['US']&",
                   "fields=id,page_id,ad_delivery_start_time,ad_delivery_stop_time,spend,region_distribution,demographic_distribution,impressions&",
                   "limit={num}") %>%
        xml2::url_escape(reserved="=&/:?,|")
}

parse_resp = function(resp) {
    ensure_cols = c(ad_delivery_stop_time=NA_character_, ad_delivery_start_time=NA_character_)
    resp$data %>%
        as_tibble %>%
        add_column(!!!ensure_cols[setdiff(names(ensure_cols), names(.))]) %>%
        transmute(id = id,
                  start_time = lubridate::ymd_hms(ad_delivery_start_time),
                  end_time = lubridate::ymd_hms(ad_delivery_stop_time),
                  candidate = names(page_ids)[match(page_id, page_ids)],
                  cost_low = as.numeric(spend$lower_bound),
                  cost_high = as.numeric(spend$upper_bound),
                  cost = sqrt((cost_low+1) * cost_high),
                  views_low = as.numeric(impressions$lower_bound),
                  views_high = as.numeric(impressions$upper_bound),
                  views = sqrt((views_low+1) * views_high),
                  regions = map(region_distribution, ~ as_tibble(.)),
                  demgs = map(demographic_distribution, ~ as_tibble(.)))
}

get_data = function(start_date=Sys.Date()-7, end_date=Sys.Date()) {
    if (is.character(start_date)) start_date = lubridate::as_date(start_date)
    if (is.character(end_date)) end_date = lubridate::as_date(end_date)
    i = 1
    cat("Page ", i, "...\n", sep=""); i = i+1;
    resp = suppressWarnings(fromJSON(get_url(page_ids, 800), simplifyVector=T))
    d = parse_resp(resp)
    while (min(d$start_time) >= start_date && !is.null(resp$paging[["next"]])) {
        cat("Page ", i, "...\n", sep=""); i = i+1;
        resp = suppressWarnings(fromJSON(resp$paging[["next"]], simplifyVector=T))
        d = bind_rows(d, parse_resp(resp))
    }
    filter(d, start_time >= start_date, coalesce(end_time <= end_date+1, T))
}

d = get_data()

by_region = d %>%
    select(id, candidate, starts_with("cost"), regions) %>%
    unnest(regions) %>%
    rename(state=region) %>%
    group_by(candidate, state) %>%
    summarize(spending = sum(cost * as.numeric(percentage))) %>%
    left_join(state_abbr, by="state") %>%
    left_join(state_pop, by=c("abbr"="state")) %>%
    left_join(states_2016, by="abbr") %>%
    mutate(spending_pc = spending / votes)


state_d = read_csv("docs/state_history.csv") %>%
    filter(date == max(date)) %>%
    select(abbr=state, ev, prob, dem_exp, tipping_pt, rel_voter_power) %>%
    right_join(by_region, by="abbr")

state_d %>%
    pivot_wider(names_from=candidate, values_from=c(spending, spending_pc)) %>%
ggplot(aes(1000*spending_pc_biden, 1000*spending_pc_trump, size=rel_voter_power, 
           color=prob, label=abbr)) +
    geom_point() +
    geom_text(size=2, color="black") +
    scale_x_continuous(labels=scales::dollar, trans="log") +
    scale_y_continuous(labels=scales::dollar, trans="log") +
    scale_color_gradient2(midpoint=0.5)

state_d %>%
    pivot_wider(names_from=candidate, values_from=c(spending, spending_pc)) %>%
ggplot(aes(tipping_pt+1e-3, spending_biden, size=rel_voter_power, 
           color=prob, label=abbr)) +
    geom_point() +
    geom_text(size=2, color="black") +
    scale_x_continuous(labels=scales::percent, trans="logit") +
    scale_y_continuous(labels=scales::dollar, trans="log", limits=c(50, NA)) +
    geom_smooth(method=lm) +
    scale_color_gradient2(midpoint=0.5)

state_d %>%
    pivot_wider(names_from=candidate, values_from=c(spending, spending_pc)) %>%
ggplot(aes(spending_biden/spending_trump, tipping_pt+1e-4, size=rel_voter_power, 
           color=prob, label=abbr)) +
    geom_point() +
    geom_text(size=2, color="black") +
    scale_y_continuous(labels=scales::percent, trans="log") +
    scale_x_continuous(trans="log") +
    geom_smooth(method=lm) +
    scale_color_gradient2(midpoint=0.5)

state_d %>%
    pivot_wider(names_from=candidate, values_from=c(spending, spending_pc)) %>%
    filter(tipping_pt > 0) %>%
    lm(log(spending_trump) ~ qlogis(tipping_pt) + qlogis(prob) + log(votes), data=.) %>%
    anova
    summary
state_d %>%
    pivot_wider(names_from=candidate, values_from=c(spending, spending_pc)) %>%
    filter(tipping_pt > 0) %>%
    lm(log(spending_biden) ~ qlogis(tipping_pt) + qlogis(prob) + log(votes), data=.) %>%
    summary

state_d %>%
    pivot_wider(names_from=candidate, values_from=c(spending, spending_pc)) %>%
    filter(rel_voter_power > 0) %>%
    lm(log(spending_pc_trump) ~ log(rel_voter_power) + qlogis(prob), data=.) %>%
    summary 
    qplot(qlogis(1-prob), log(spending_pc_trump), data=.)
    
state_d %>%
    pivot_wider(names_from=candidate, values_from=c(spending, spending_pc)) %>%
    qplot(qlogis(1e-5+prob), log(spending_pc_trump), data=.)
