get_natl_prior_m = function(path, recompile=F) {
    compiled_model = paste0(path, ".rdata")
    if (recompile && file.exists(compiled_model))
        file.remove(compiled_model)

    if (file.exists(compiled_model)) {
        model = read_rds(compiled_model)
    } else {
        natl_model_d = suppressWarnings(suppressMessages(
            read_csv("data/historical/tfc.csv")
            )) %>%
            transmute(year = year, 
                      q2_gdp = q2_growth/100, 
                      q2_appr = qlogis(0.5 + approval/200),
                      two_term = incumbent,
                      inc_vote = qlogis(vote/100))
        model = stan_lm(inc_vote ~ q2_gdp + q2_appr + two_term, data=natl_model_d, 
                        prior=R2(location=0.8), chains=1, iter=2000, 
                        control=list(adapt_delta=0.999))
        
        write_rds(model, compiled_model, compress="gz")
    }

    model
}

get_state_prior_m = function(path, state_d, recompile=F) {
    compiled_model = paste0(path, ".rdata")
    if (recompile && file.exists(compiled_model))
        file.remove(compiled_model)

    if (file.exists(compiled_model)) {
        model = read_rds(compiled_model)
    } else {
        d = state_d %>%
            group_by(state) %>%
            mutate(dem = qlogis(dem / (dem + gop)) - qlogis(natl_dem),
                   lpr_state = lag(pr_state, order_by=year),
                   l2pr_state = lag(pr_state, 2, order_by=year),
                   lvp_state = lag(vp_state, order_by=year),
                   ldem = lag(dem, order_by=year),
                   l2dem = lag(dem, 2, order_by=year),
                   inc = 2*(inc.party=="D") - 1,
                   linc = lag(inc, order_by=year),
                   ltt = lag(two.term, order_by=year),
                   regn_yr = str_c(region, "_", year %% 100),
                   q2.appr = qlogis(q2.appr)) %>%
            drop_na
        
        model = stan_lmer(dem ~ ldem + l2dem + (ldem+l2dem|region) + 
                              (1|region:year) + inc*ldem*two.term + 
                              pr_state + vp_state + lpr_state + l2pr_state + lvp_state,
                          data=d, chains=1, warmup=500, iter=1700, prior=cauchy())
        
        write_rds(model, compiled_model, compress="gz")
    }

    model
}


get_polls_m = function(path, recompile=F) {
    compiled_model = path
    if (recompile && file.exists(compiled_model))
        file.remove(compiled_model)

    cmdstan_model("stan/polls.stan")#, compiler_flags="CXXFLAGS += -Ofast")
}

