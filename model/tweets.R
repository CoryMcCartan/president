suppressMessages(library(tidyverse))
suppressMessages(library(rpredictit))
suppressMessages(library(lubridate))
suppressMessages(library(slider))
suppressMessages(library(BART))
suppressMessages(library(rstanarm))

raw_tweets = read_rds("data/tweets.rdata") %>%
    select(tstamp=created_at, rt=is_retweet) %>%
    mutate(tstamp = with_tz(mdy_hms(tstamp, tz="UTC"), "America/New_York"),
           hour = (hour(tstamp) + 12) %% 24,
           date = as_date(tstamp),
           date = if_else(hour >= 12, date - 1, date),
           date = force_tz(as_datetime(date+1),
                           "America/New_York") - 12*3600,
           rt = coalesce(rt, F))

start_time = as_date(min(raw_tweets$date))
end_time = as_date(max(raw_tweets$date))
dates = force_tz(as_datetime(seq(start_time, end_time, by="day")+1),
                 tzone="America/New_York") - 12*3600

tweets = raw_tweets %>%
    group_by(date, hour) %>%
    summarize(n = n(),
              rt = mean(rt)) %>%
    ungroup
missing = crossing(hour=0:23, date=dates) %>%
    anti_join(tweets, by=c("date", "hour")) %>%
    arrange(date, hour)
tweets = full_join(tweets, missing, by=c("date", "hour")) %>%
    arrange(date, hour) %>%
    mutate(n = coalesce(n, 0L), rt = coalesce(rt, 0),
           week = floor_date(date, "week", week_start=3) + 12*3600,
           day = wday(date, week_start=3)) %>%
    select(week, date, day, hour, n, rt) %>%
    group_by(week) %>%
    filter(head(day, 1) == 1, tail(day, 1) == 7) %>%
    mutate(n_cuml = cumsum(n),
           n_rem = sum(n) - n_cuml) %>%
    ungroup
pattern = tweets %>%
    select(week, day, hour, n, n_cuml) %>%
    group_by(week) %>%
    arrange(week, day, hour) %>%
    mutate(total = sum(n),
        frac = (n_cuml - n/2)/ sum(n)) %>%
    group_by(day, hour) %>%
    summarize(mean_rem = 1-mean(frac))

tweets = tweets %>%
    left_join(pattern, by=c("day", "hour")) %>%
    mutate(lpace = log(n_cuml+1) - log(1 - mean_rem),
           lpace_fwd = log(n_rem+1) - log(mean_rem),
           n_h1 = slide_vec(n, sum, .before=1),
           n_h6 = slide_vec(n, sum, .before=6),
           n_h12 = slide_vec(n, sum, .before=12),
           n_d1 = slide_vec(n, sum, .before=1*24),
           n_d2 = slide_vec(n, sum, .before=2*24),
           n_d3 = slide_vec(n, sum, .before=3*24),
           n_d4 = slide_vec(n, sum, .before=4*24),
           n_d5 = slide_vec(n, sum, .before=5*24),
           n_d6 = slide_vec(n, sum, .before=6*24),
           n_d7 = slide_vec(n, sum, .before=7*24),
           n_d10 =slide_vec(n, sum, .before=10*24),
           n_w2 = slide_vec(n, sum, .before=2*7*24),
           n_w3 = slide_vec(n, sum, .before=3*7*24),
           n_w4 = slide_vec(n, sum, .before=4*7*24),
           n_w6 = slide_vec(n, sum, .before=6*7*24),
           n_w8 = slide_vec(n, sum, .before=8*7*24),
           rt_d1 = slide_vec(n*rt, sum, .before=1*24) / n_d1,
           rt_d7 = slide_vec(n*rt, sum, .before=7*24) / n_d7,
           rt_w8 = slide_vec(n*rt, sum, .before=8*7*24) / n_w8) %>%
    filter(date >= start_time + 8*7 + 1)


ggplot(tweets, aes(lpace, lpace_fwd, color=24*day+hour)) +
    geom_point(alpha=0.2, size=0.2)

ggplot(tweets, aes(date+hour*3600, n)) + geom_line()
ggplot(tweets, aes(date+hour*3600, n_d10)) + geom_line()
ggplot(tweets, aes(date+hour*3600, n_w8)) + geom_line()
ggplot(tweets, aes(date+hour*3600, n_rem)) + geom_line()

test_date = ymd("2019-09-01")
d_train = filter(drop_na(tweets), week <= test_date)
d_test = filter(drop_na(tweets), week > test_date)

m = lm(lpace_fwd ~ log(n_cuml-n+1) + log(n_h1+1) + log(n_h6+1) +
            log(n_h12+1) + log(n_d1+1) + log(n_d2+1) + log(n_d3) +
            log(n_d5) +  log(n_d7) + log(n_d10) + log(n_w2) + log(n_w3) +
            log(n_w4) + log(n_w6) + log(n_w8) + rt_d7 + rt_w8 + I(hour^2) +
            as.factor(day)*hour + lpace,
       data=d_train)
d_test$y_pred = exp(predict(m, newdata=d_test) + log(d_test$mean_rem)) - 1
d_test$y_pred = predict(m, newdata=d_test, type="response")
with(d_test, sd(y_pred - n_rem))
with(d_test, median(abs(y_pred - n_rem)))
with(d_test, cor(y_pred, n_rem)^2)

d_test %>%
    mutate(error = y_pred - n_rem) %>%
    group_by(day, hour) %>%
    summarize(se=sd(error), mad=mean(abs(error)),
              rsq=100*cor(n_rem, y_pred)^2,
              rsq_ref = 100*cor(n_cuml, n_rem)^2) %>%
    pivot_longer(-day:-hour) %>%
ggplot(aes(day+hour/24, value, color=name)) + geom_line(size=1)

form = ~ I(n_cuml-n) + mean_rem + n_h1 + n_h6 + n_h12 + n_d1 + n_d2 + n_d3 +
    n_d4 + n_d5 + n_d6 + n_d7 + n_d10 + n_w2 + n_w3 + n_w4 + n_w6 + n_w8 +
    rt_d1 + rt_d7 + rt_w8 + as.factor(day) + lpace
idxs = 1:nrow(d_train)#sample.int(nrow(d_train), 5000)
x_train = model.matrix(form, data=d_train)[idxs,]
y_train = with(d_train[idxs,], lpace_fwd)
x_test = model.matrix(form, data=d_test)

mb = wbart(x_train, y_train, x_test, sparse=T, a=0.5, ntree=200)
d_test$y_pred = exp(mb$yhat.test.mean + log(d_test$mean_rem)) - 1

