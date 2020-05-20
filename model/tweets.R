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
    mutate(n_cuml = cumsum(n),
           n_rem = sum(n) - n_cuml) %>%
    ungroup %>%
    mutate(hr_rem = 7*24 - ((day - 1)*24 + hour),
           pace = hr_rem * (n_cuml + 0.5) / ((day - 1)*24 + hour + 0.5),
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
    filter(date >= start_time + 8*7)

ggplot(tweets, aes(date+hour*3600, n)) + geom_line()
ggplot(tweets, aes(date+hour*3600, n_d10)) + geom_line()
ggplot(tweets, aes(date+hour*3600, n_w8)) + geom_line()
ggplot(tweets, aes(date+hour*3600, n_rem)) + geom_line()

test_date = ymd("2019-11-01")
d_train = filter(drop_na(tweets), week <= test_date)
d_test = filter(drop_na(tweets), week > test_date)

r2_d = function(d) {
    m = stan_glm(n_rem ~ log(n_cuml+0.5)*hour + log(n_h1+0.5) + log(n_h6+0.5) +
                log(n_h12+0.5) + log(n_d1+0.5) + log(n_d2+0.5) + log(n_d3) +
                log(n_d5) +  log(n_d7) + log(n_d10) + log(n_w2) + log(n_w3) +
                log(n_w4) + log(n_w6) + log(n_w8) + rt_d7 + rt_w8 + I(hour^2),
            data=filter(d_train, day==d),
            prior=cauchy(), #chains=1,
            algorithm="meanfield",
            #family=quasipoisson())
            family=neg_binomial_2())

    d_test_d = filter(d_test, day==d)
    y_pred = suppressWarnings(predict(m, newdata=d_test_d, type="response"))
    y_act = d_test_d$n_rem
    sd(y_pred - y_act)
    cor(y_pred, y_act)^2
    median(abs(y_pred - y_act))
}
map_dbl(1:7, r2_d)

m = lm(log(n_rem+1) ~ log(n_cuml+1) + log(n_h1+1) + log(n_h6+1) +
            log(n_h12+1) + log(n_d1+1) + log(n_d2+1) + log(n_d3) +
            log(n_d5) +  log(n_d7) + log(n_d10) + log(n_w2) + log(n_w3) +
            log(n_w4) + log(n_w6) + log(n_w8) + rt_d7 + rt_w8 + I(hour^2) +
            as.factor(day)*hour + log(pace),# + (1|week),
       data=d_train)#,
       #prior=cauchy(), chains=1)
#y_pred = colMeans(posterior_predict(m, newdata=d_test))
d_test$y_pred = exp(predict(m, newdata=d_test, allow.new.levels=T)) - 1
d_test$y_pred = predict(m, newdata=d_test, type="response")
with(d_test, sd(y_pred - n_rem))
with(d_test, cor(y_pred, n_rem)^2)
with(d_test, median(abs(y_pred - n_rem)))

d_test %>%
    mutate(error = y_pred - n_rem) %>%
    group_by(day, hour) %>%
    summarize(se=sd(error), mad=mean(abs(error)),
              rsq=100*cor(n_rem, y_pred)^2,
              rsq_ref = 100*cor(n_cuml, n_rem)^2) %>%
    pivot_longer(-day:-hour) %>%
ggplot(aes(day+hour/24, value, color=name)) + geom_line(size=1)


select(tail(d_test, 20), date, day, hour, n, n_cuml, n_rem, y_pred, hr_rem, pace) %>%
    mutate(date = date + hour*3600,
           wk  = n_cuml + y_pred)

d_test %>%
    mutate(error = y_pred - n_rem) %>%
ggplot(aes(error)) +
    facet_wrap(~ day) +
    geom_histogram(bins=10)

m = glm(n_rem ~ n_cuml + n_h1 + n_h6 +
            n_h12 + n_d1 + n_d2 + n_d3 +
            n_d5 +  n_d7 + n_d10 + n_w2 + n_w3 +
            n_w4 + n_w6 + n_w8 + rt_d7 + rt_w8 + I(hour^2) +
            as.factor(day)*hour + pace,
       data=d_train,
       family=quasipoisson())



x_train = model.matrix(form, data=d_train)
y_train = d_train$n_rem
x_test = model.matrix(form, data=d_test)
y_test = d_test$n_rem
form = n_rem ~ n_cuml + n_h1 + n_h6 + n_h12 + n_d1 + n_d2 + n_d3 + n_d4 + n_d5 +
    n_d6 + n_d7 + n_d10 + n_w2 + n_w3 + n_w4 + n_w6 + n_w8 + rt_d1 + rt_d7 +
    rt_w8 + as.factor(day)

mb = wbart(x_train, y_train, x_test, sparse=T, a=0.6, ntree=150)

cor(mb$yhat.test.mean, y_test)^2
sd(mb$yhat.test.mean - y_test)
median(abs(mb$yhat.test.mean - y_test))

n_rem ~ n_cuml + n_h1 + n_h6 + n_h12 + n_d1 + n_d2 + n_d3 + n_d4 + n_d5 +
    n_d6 + n_d7 + n_d10 + n_w2 + n_w3 + n_w4 + n_w6 + n_w8 + rt_d1 + rt_d7 +
    rt_w8 + as.factor(day)
