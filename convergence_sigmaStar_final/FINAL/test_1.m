% Experimental result for crosses and circles (1+1)-ES
% Use simulation


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

i = 1;
for k = 1:1:V_LENGTH 
    v_temp = v_array(k); 
    sigma_star = opt_sigma_star(k);
    
    sigma_ep_star = sigma_ep_star_array(k);
    for j = 1:1:NUM_OF_RUNS
        x0 = randn(10,1); 
        a = one_plus_one_noise(f,x0,sigma_star,sigma_ep_star,NUM_OF_ITERATIONS);
        
        x0 = randn(100,1); 
        b = one_plus_one_noise(f,x0,sigma_star,sigma_ep_star,NUM_OF_ITERATIONS);

        temp_sigma_f_x_array_10(j,:) = cell2mat(a(3));
        temp_sigma_f_x_array_100(j,:) = cell2mat(b(3));
        % avoid nan in fx (if nan appears -> set 0)
        if(sum(isnan(temp_sigma_f_x_array_10(j,:))) == 0)
            temp_sigma_counvergence_rate_array_10(j) = cell2mat(a(7));
        else 
            temp_sigma_counvergence_rate_array_10(j) = 0;            
        end
        if(sum(isnan(temp_sigma_f_x_array_100(j,:))) == 0)
            temp_sigma_counvergence_rate_array_100(j) = cell2mat(b(7));
        else
            temp_sigma_counvergence_rate_array_100(j) = 0;
        end
    end
    % STORE averaged MAX convergence rate 
    sigma_counvergence_rate_array_10(i) = median(temp_sigma_counvergence_rate_array_10);
    sigma_counvergence_rate_array_100(i) = median(temp_sigma_counvergence_rate_array_100);
    
    if(k == V_LENGTH)
        disp(sigma_counvergence_rate_array_100);
    end
    i = i + 1;
        
end

scatterColour = 'r';
FIG_NUM = 1;
subplotNum_step = 2;

hold on
legend('-DynamicLegend'); 

figure(FIG_NUM);hold on;
subplot(1,2,subplotNum_step-1);hold on;
n=10;
typeDot = 'x';
d =sprintf('N = %.0f',n);
scatter(v_array, sigma_counvergence_rate_array_10,typeDot,scatterColour,'DisplayName',d); hold on; 
legend('-DynamicLegend'); 
legend('show');
leg = legend();
leg.FontSize = 10;
title(leg,'Dimension of data');

subplot(1,2,subplotNum_step-1);hold on;
n=100;
typeDot = 'o';
d =sprintf('N = %.0f',n);
scatter(v_array, sigma_counvergence_rate_array_100,typeDot,scatterColour,'DisplayName',d); hold on; 
legend('-DynamicLegend'); 
legend('show');
xlabel('noise-to-signal ratio \vartheta','fontsize',20);
set(gca, 'XScale', 'log')
ylabel('opt. expected fitness gain \eta_{opt}','FontSize',15); 
xlim([0.1 10]);
ylim([0 2]);

