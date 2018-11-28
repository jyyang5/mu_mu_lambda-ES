f = @(x) (x'*x);        % quadratic sphere

NUM_OF_RUNS = 50;
NUM_OF_ITERATIONS = 10000;

v_expedted_curve_array = exp(-2.302585092994046: 0.0461:2.302585092994046+0.01);
v_array = v_expedted_curve_array(1:10:101);
V_LENGTH = length(v_array);
opt_sigma_star = [7.6100, 5.8100, 4.3700,3.2500,2.4400,1.9000,1.5900,1.4200,1.3300,1.2900,1.2600];

sigma_ep_star_array = v_array.*opt_sigma_star;

eta_10 = zeros(11,1);
eta_100 = zeros(11,1);

sigma_counvergence_rate_array_10 = zeros(1,V_LENGTH);
sigma_counvergence_rate_array_100 = zeros(1,V_LENGTH);

k = 1;
sigma_star = opt_sigma_star(k);
sigma_ep_star = sigma_ep_star_array(k);

x0 = randn(100,1); 
b = one_plus_one_noise(f,x0,sigma_star,sigma_ep_star,NUM_OF_ITERATIONS);