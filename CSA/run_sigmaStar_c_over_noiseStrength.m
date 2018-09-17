
% run sigma* and convergence rate over sigma_ep_star
f = @(x) (x'*x);

mu = 3;
n = 10;                 % dim of the data
lambda = 10;
x0 = randn(n,mu); 

NUM_OF_ITERATIONS = 4000;
NUM_OF_RUNS = 40;




% matrix for different sigma_ep_star matrix(sigma_i,median of NUM_OF_RUNS)
s_start = 0.2;
increment =0.2;
s_end = 8;


temp_sigma_sigma_star_array = zeros(1,NUM_OF_RUNS);
temp_sigma_T_array = zeros(1,NUM_OF_RUNS);
temp_sigma_f_x_array = zeros(NUM_OF_RUNS,10000);
temp_sigma_counvergence_rate_array = zeros(1,NUM_OF_RUNS);

sigma_sigma_star_array = zeros(1,(s_end-s_start)/increment+1);
sigma_T_array = zeros(1,(s_end-s_start)/increment+1);
sigma_f_x_array = zeros((s_end-s_start)/increment+1,10000);
sigma_counvergence_rate_array = zeros(1,(s_end-s_start)/increment+1);

i = 1;
for sigma_ep_star = s_start:increment:s_end
    for j = 1:1:NUM_OF_RUNS
        a = mml_csa(f,x0,sigma0,sigma_ep_star,lambda,NUM_OF_ITERATIONS);
        temp_sigma_f_x_array(j,:) = cell2mat(a(3));
        temp_sigma_sigma_star_array(j) = cell2mat(a(8));
        temp_sigma_counvergence_rate_array(j) = cell2mat(a(7));
        temp_sigma_T_array(j) = cell2mat(a(1));

    end
    for j = 1:1:10000
        sigma_f_x_array(i,j) = median(temp_sigma_f_x_array(:,j));
    end
    sigma_sigma_star_array(i) = median(temp_sigma_sigma_star_array);
    sigma_counvergence_rate_array(i) = median(temp_sigma_counvergence_rate_array);
    sigma_T_array(i) = median(temp_sigma_T_array);
    
    
    i = i + 1;
    
end

figure(1)
subplot(1,2,1)
plot(s_start:increment:s_end, sigma_counvergence_rate_array);
ylabel('convergnece rate c','fontsize',20);
xlabel('noise strength \sigma_\epsilon*','fontsize',20);
%set(gca,'xscale','log')

subplot(1,2,2)
plot(s_start:increment:s_end, sigma_sigma_star_array);
ylabel('normalzied step size \sigma*','fontsize',20);
xlabel('noise strength \sigma_\epsilon*','fontsize',20);
