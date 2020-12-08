data{
    int W[1345];
    int S[1345];
    int B[1345];
    int C[1345];
    int L[1345];
    int I[1345];
}
parameters{
    vector[17] a;
    vector[4] b;
    vector[5] g;
    vector[2] d;
    vector[3] e;
    real a_bar;
    real<lower=0> sigma_a;
    real<lower=0> sigma_b;
    real<lower=0> sigma_g;
    real<lower=0> sigma_d;
    real<lower=0> sigma_e;
}
model{
    vector[1345] p;
    sigma_e ~ exponential( 1 );
    sigma_d ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_a ~ exponential( 1 );
    a_bar ~ normal( 0 , 1.5 );
    e ~ normal( 0 , 1 );
    d ~ normal( 0 , 1 );
    g ~ normal( 0 , 1 );
    b ~ normal( 0 , 1 );
    a ~ normal( 0 , 1 );
    for ( i in 1:1345 ) {
        p[i] = a_bar + a[I[i]] * sigma_a + b[L[i]] * sigma_b + g[C[i]] * sigma_g + d[B[i]] * sigma_d + e[S[i]] * sigma_e;
        p[i] = inv_logit(p[i]);
    }
    W ~ binomial( 1 , p );
}
generated quantities{
    vector[1345] log_lik;
    vector[1345] p;
    for ( i in 1:1345 ) {
        p[i] = a_bar + a[I[i]] * sigma_a + b[L[i]] * sigma_b + g[C[i]] * sigma_g + d[B[i]] * sigma_d + e[S[i]] * sigma_e;
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:1345 ) log_lik[i] = binomial_lpmf( W[i] | 1 , p[i] );
}


