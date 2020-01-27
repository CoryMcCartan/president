library(tidyverse)
library(lubridate)
library(tidybayes)

election = ymd("2016-11-08")
start_date = election - 20*7

d = suppressMessages(read_csv("data/polls/president_general_polls_2016.csv")) %>%
    filter(type == "polls-only") %>%
    mutate(state = na_if(state, "U.S."),
           national = is.na(state),
           date = mdy(startdate) + (mdy(enddate) - mdy(startdate))/2,
           day = as.numeric(date - start_date),
           week = as.numeric(floor(day / 7)),
           dem = rawpoll_clinton / (rawpoll_clinton + rawpoll_trump),
           gop = rawpoll_trump / (rawpoll_clinton + rawpoll_trump),
           pollster = abbreviate(pollster, 8, named=F),
           sample = if_else(population=="v", "rv", population)) %>%
    filter(day >= 0) %>%
    select(state, national, day, week, pollster, size=samplesize, sample, dem, gop)

model.d = compose_data(d, 
                       D=length(unique(d$day)), 
                       W=length(unique(d$week)), 
                       )
