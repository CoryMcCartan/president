data {
    int N;
    int T;
    int S;
    int K;
    int<lower=0, upper=1> estimate;
    
    vector[N] inc_party;
    int year[N];
    int state[N];
    matrix[N, K] X;
    vector[N] y;
    
    real prior_theta_scale;
    real prior_sigma_e_scale;
    real prior_sigma_t_scale;
    real prior_sigma_s_scale;
    real prior_sigma_s_location;
    real prior_eta_scale;
}
transformed data {
    // QR decomp. of covariates
    matrix[N, K] Q_ast;
    matrix[K, K] R_ast;
    matrix[K, K] R_ast_inverse;
    //vector[N] year_v;
    //vector[N] year_v_sq;
    
    Q_ast = qr_thin_Q(X);
    R_ast = qr_thin_R(X);
    R_ast_inverse = inverse(R_ast);
    
    //year_v = to_vector(year);
    //year_v_sq = square(year_v);
}
parameters {
    //vector[T] a_t_std;
    vector[T] a_t;
    vector[K] theta;
    //matrix[T,S] a_s_std;
    matrix[T,S] a_s;
    //vector[S] a_s;
    //vector[S] beta_s;
    //vector[S] beta_2s;
    
    real<lower=0> sigma_e;
    real<lower=0> sigma_t;
    real<lower=0> sigma_s;
    //real<lower=0,upper=0.05> eta_s;
}
transformed parameters {
    vector[N] lin_pred;
    //matrix[T,S] a_s;
    //vector[T] a_t;
    
    //a_t = a_t_std * sigma_t;
    
    // investigate whether dropping sigma_s made a difference
    //a_s[1] = a_s_std[1] * sigma_s;
    //for (t in 2:T)
    //    a_s[t] = a_s[t-1] + a_s_std[t]*prior_eta_scale;
            
    for (j in 1:N) 
        lin_pred[j] = a_t[year[j]] + Q_ast[j]*theta 
            + a_s[year[j],state[j]]*inc_party[j];
            
    //y_pred = a_t[year] + Q_ast*theta  + (a_s[state] 
    //    + beta_s[state].*year_v + beta_2s[state].*year_v_sq).*inc_party;
}
model {
    if (estimate == 1)
        logit(y) ~ normal(lin_pred, sigma_e);
    
    //a_t_std ~ std_normal();        
    //for (t in 1:T)
    //    a_s_std[t] ~ std_normal();
    a_t ~ normal(0, sigma_t);        
    a_s[1] ~ normal(0, sigma_s);
    for (t in 2:T)
        a_s[t] ~ normal(a_s[t-1], prior_eta_scale);
    //a_s ~ normal(0, sigma_s);
    
    theta ~ normal(0, prior_theta_scale);
    sigma_e ~ normal(0, prior_sigma_e_scale);
    sigma_t ~ normal(0, prior_sigma_t_scale);
    sigma_s ~ normal(prior_sigma_s_location, prior_sigma_s_scale);
    //beta_s ~ normal(0, prior_eta_scale); 
    //beta_2s ~ normal(0, prior_eta_scale); 
    //eta_s ~ normal(0, prior_eta_scale);
}
generated quantities {
    vector[K] beta;
    vector[N] log_lik;
    vector[N] y_pred;
    
    beta = R_ast_inverse * theta;
    for (j in 1:N) {
        log_lik[j] = normal_lpdf(logit(y[j]) | lin_pred[j], sigma_e);
        y_pred[j] = inv_logit(normal_rng(lin_pred[j], sigma_e));
    }
}
