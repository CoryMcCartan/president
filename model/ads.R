library(tidyverse)
library(jsonlite)
library(usmap)

state_abbr = read_rds("output/state_data_2016.rdata") %>%
    select(state, abbr) %>%
    filter(!str_detect(state, "CD-"))
state_abbr$state[state_abbr$abbr=="DC"] = "Washington, District of Columbia"

state_pop = suppressMessages(read_csv("data/state-returns/1976-2016-president.csv")) %>% 
    filter(year==2016) %>% 
    select(state=state_po, votes=totalvotes) %>%
    distinct()

page_ids = tibble(`Page ID`=c(7860876103, 153080620724), candidate=c("Biden", "Trump"))

{
pb = progress_estimated(51)
d = map_dfr(state_abbr$state, function(s) {
    pb$tick()$print()
    url = glue("data/ads/report/regions/FacebookAdLibraryReport_",
            "2020-03-31_US_last_30_days_{s}.csv")
    suppressMessages(read_csv(url)) %>%
        inner_join(page_ids, by="Page ID") %>%
        select(candidate, spending=4) %>%
        group_by(candidate) %>%
        summarize(state = s,
                  spending = sum(as.numeric(spending)))
}) %>%
    left_join(state_abbr, by="state") %>%
    left_join(state_pop, by=c("abbr"="state")) %>%
    select(state, abbr, votes, candidate, spending) %>%
    mutate(spending_pc = spending / votes)
d$state[d$abbr == "DC"] = "District of Columbia"
}

#plot_usmap("states", value="spending_pc", data=filter(d, candidate=="Trump")) + guides(fill=F)


state_d = read_csv("docs/state_history.csv") %>%
    filter(date == max(date)) %>%
    select(abbr=state, ev, prob, dem_exp, tipping_pt, rel_voter_power) %>%
    right_join(d, by="abbr")

state_d %>%
    mutate(candidate = str_to_lower(candidate)) %>%
    pivot_wider(names_from=candidate, values_from=c(spending, spending_pc)) %>%
ggplot(aes(spending_pc_biden, spending_pc_trump, size=rel_voter_power, 
           color=prob, label=abbr)) +
    geom_point() +
    geom_text(size=2, color="black") +
    scale_x_continuous(labels=scales::dollar, trans="log") +
    scale_y_continuous(labels=scales::dollar, trans="log") +
    scale_color_gradient2(midpoint=0.5)

state_d %>%
    mutate(label = paste(abbr, candidate)) %>%
    filter(tipping_pt > 0.01) %>%
ggplot(aes(tipping_pt, spending, group=abbr, color=candidate, label=abbr)) +
    geom_line(color="#00000066") +
    geom_text(size=2, fontface="bold") +
    scale_color_manual(values=c("#6677ff", "#ff7766")) +
    scale_y_continuous(labels=scales::dollar, trans="log")

    