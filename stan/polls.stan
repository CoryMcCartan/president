/************************************************
 * PRESIDENTIAL ELECTION POLLING MODEL          *
 * CORY McCARTAN                                *
 * (c) 2020                                     *
 ************************************************/

data {
    int N; // number of polls
    int D; // "days"
    int W; // "weeks"
    int D_W; // "days" per "week"
    int N_state;
    int N_regn;
    int N_firm;
    
    int<lower=1> day[N];
    int<lower=1> week[N];
    int<lower=-1> state[N];
    int<lower=1> firm[N];
    int<lower=1> week_day[D];
    int<lower=1> regn[N];
    real<lower=0, upper=1> week_frac[D];
    vector<lower=0, upper=1>[N] type_rv;
    vector<lower=0, upper=1>[N] type_lv;
    vector<lower=0, upper=1>[N] type_a;
    vector<lower=0, upper=1>[N] national;
    
    vector[N] dem;
    vector[N] var_poll; 
    
    real prior_natl_mean;
    vector[N_state] prior_state_mean;
    cov_matrix[N_state] prior_state_cov; // prior model result correlation
    
    real<lower=0> prior_natl_sd;
    real prior_natl_poll_bias;
    real<lower=0> prior_natl_poll_error;
    real prior_all_state_poll_bias;
    real<lower=0> prior_all_state_poll_error; 
    real prior_regn_poll_bias;
    real<lower=0> prior_regn_poll_error;
    real<lower=0> prior_states_poll_error; 
    // sds state poll errors around global state+natl error
    
    real prior_rv_bias;
    real prior_lv_bias;
    real prior_a_bias;
    real<lower=0> lv_rv_ratio;
}

transformed data {
    vector[N] logit_dem = logit(dem);
    cholesky_factor_cov[N_state] chol_state_cov 
        = cholesky_decompose(prior_state_cov);
}

parameters {
    real<lower=0> sigma_natl; // day-to-day variance
    real<lower=0> sigma_state; // week-to-week variance
    
    real<lower=0> nonsamp_var;
    real natl_error;
    real all_state_error;
    vector[N_regn] regn_error;
    vector[N_state] state_error;
    real<lower=0> sd_firm; // hyperparameter for pollster errors
    real bias_rv; // registered voter poll bias
    real bias_lv; // likely voter poll bias
    real bias_a; // adult voter poll bias
    
    // reparametrization
    vector[D] delta_natl; // steps of random walk
    vector[N_state] delta_state[W]; // steps of randowm walk
    vector[N_firm] delta_firm; // polling firm standardized errors
}

transformed parameters {
    vector[D] mu_natl; // national voter intention, logit
    vector[N_state] mu_state[W]; // national voter intention, logit
    
    mu_natl[D] = prior_natl_mean + prior_natl_sd*delta_natl[D];
    mu_state[W] = prior_state_mean + chol_state_cov*delta_state[W];
    for (i in 1:(D-1)) {
        int d = D - i;
        mu_natl[d] = mu_natl[d+1] + sigma_natl*delta_natl[d];
    }
    for (i in 1:(W-1)) {
        int w = W - i;
        mu_state[w] = mu_state[w+1] + sqrt(D_W)*sigma_state*delta_state[w];
    }
    
}

model {
    // support for dem. in specific poll
    vector[N] val_poll = mu_natl[day] + 4*sd_firm*delta_firm[firm] 
        + 4*bias_rv*type_rv + 4*bias_lv*type_lv + 4*bias_a*type_a + natl_error;
    for (i in 1:N) {
        if (!national[i]) {
            val_poll[i] += all_state_error + regn_error[regn[i]] 
                + prior_states_poll_error*state_error[state[i]];
            if (day[i] <= D_W)
                val_poll[i] += mu_state[week[i], state[i]];
            else
                val_poll[i] += week_frac[day[i]]*mu_state[week[i], state[i]]
                    + (1 - week_frac[day[i]])*mu_state[week[i] - 1, state[i]];
        }
    }
    
    logit_dem ~ normal(val_poll, sqrt(var_poll + nonsamp_var)); // modified binom. approx.
    
    // reparametrization
    delta_natl ~ std_normal();
    for (w in 1:W) {
        delta_state[w] ~ std_normal();
    }
    delta_firm ~ std_normal();
    
    // priors
    sigma_natl ~ gamma(3, 3/0.06);
    sigma_state ~ gamma(2, 2/0.02);
    natl_error ~ normal(prior_natl_poll_bias, prior_natl_poll_error);
    all_state_error ~ normal(prior_all_state_poll_bias, prior_all_state_poll_error);
    regn_error ~ normal(prior_regn_poll_bias, prior_regn_poll_error);
    state_error ~ std_normal();
    
    sd_firm ~ gamma(2, 2/0.05);
    bias_rv ~ normal(prior_rv_bias, 0.04);
    bias_lv ~ normal(prior_lv_bias, 0.04/lv_rv_ratio);
    bias_a ~ normal(prior_a_bias, 0.04);
}

generated quantities {
    real nonsamp_sd = sqrt(nonsamp_var);
    vector[D] natl_dem = inv_logit(mu_natl);
    vector[N_firm] house_effects = sd_firm * delta_firm;
    vector[N_state] state_dem[D];
    
    for (d in 1:D_W) {
        state_dem[d] = inv_logit(mu_natl[d] + mu_state[week_day[d]]);
    }
    for (d in (D_W+1):D) {
        state_dem[d] = inv_logit(mu_natl[d] + week_frac[d]*mu_state[week_day[d]]
            + (1 - week_frac[d])*mu_state[week_day[d] - 1]);
    }
}
