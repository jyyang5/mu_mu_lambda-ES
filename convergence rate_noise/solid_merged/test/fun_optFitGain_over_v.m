% plot opt. normalized step size & opt. expected fitness gain over v 

function fun_optFitGain_over_v(NUM_OF_RUNS,mu,lambda,v_array,n,scatterColour,typeDot,FIG_NUM,c_mu_lambda,v_curve_array)
%Input
%   v_array:              array noise-to-signal ratio = sigma_ep_star /sigma_star
%   n:                    dim of data 
%   scatterColour:        colour of the plots
%   typeDot:              type of dot (+  Plus sign, o  Circle, *  Asterisk, x    Cross)
%   LINE_OR_NOT:          if just scatter = 0
%   c_mu_lambda:          a factor for computing n-> infty convergence rate
%   sigma_s_i_e:          [sigma_start, increment, sigma_end]
%Return
%   scatter (experimental result) or curve (expected)


% v_array = [0.1 0.25 0.4 1 2 4 8 16];
V_LENGTH = length(v_array);



f = @(x) (x'*x);

%n = 10;                 % dim of the data
%ita = 4;     
x0 = randn(n,mu); 

NUM_OF_ITERATIONS = 2000;
% NUM_OF_RUNS = 1;



% matrix for different sigma* matrix(sigma_i,median of NUM_OF_RUNS)
s_start = 0.2;
increment =0.2;
s_end = 32.2;

temp_sigma_success_rate_array = zeros(1,NUM_OF_RUNS);
temp_parent_not_lowest_of_quadratic_array = zeros(1,NUM_OF_RUNS);
temp_sigma_T_array = zeros(1,NUM_OF_RUNS);
temp_sigma_f_x_array = zeros(NUM_OF_RUNS,6000);
temp_sigma_counvergence_rate_array = zeros(1,NUM_OF_RUNS);

sigma_success_rate_array = zeros(1,(s_end-s_start)/increment+1);
parent_not_lowest_of_quadratic_array = zeros(1,(s_end-s_start)/increment+1);
sigma_T_array = zeros(1,(s_end-s_start)/increment+1);
sigma_f_x_array = zeros((s_end-s_start)/increment+1,6000);
sigma_counvergence_rate_array = zeros(1,(s_end-s_start)/increment+1);

v_convergence_rate_array = zeros(1,V_LENGTH);                               % store the max median c over v

for k = 1:1:V_LENGTH 
    v_temp = v_array(k); 
    i = 1;
    for sigma_star = s_start:increment:s_end
        for j = 1:1:NUM_OF_RUNS
            sigma_ep_star = v_temp*sigma_star;
            a = mml_noise(f,x0,sigma_star,sigma_ep_star,lambda,NUM_OF_ITERATIONS);
            %temp_parent_not_lowest_of_quadratic_array(j) = cell2mat(a(6));
            %temp_sigma_counvergence_rate_array(j) = cell2mat(a(7));
            %temp_sigma_T_array(j) = cell2mat(a(1));
%             temp_sigma_success_rate_array(j) = cell2mat(a(5));
            temp_sigma_f_x_array(j,:) = cell2mat(a(3));
            % avoid nan in fx (if nan appears -> set 0)
            if(sum(isnan(temp_sigma_f_x_array(j,:))) == 0)
                temp_sigma_counvergence_rate_array(j) = cell2mat(a(7));
                temp_sigma_T_array(j) = cell2mat(a(1));
                %temp_sigma_success_rate_array(j) = cell2mat(a(5));
            else 
                temp_sigma_counvergence_rate_array(j) = 0;
                temp_sigma_T_array(j) = 0;
                %temp_sigma_success_rate_array(j) = 0;
            end
        end
        sigma_counvergence_rate_array(i) = median(temp_sigma_counvergence_rate_array);
        sigma_T_array(i) = median(temp_sigma_T_array);
        i = i + 1;
        
    end
    v_convergence_rate_array(k) = max(sigma_counvergence_rate_array);
end

figure(FIG_NUM);
hold on
legend('-DynamicLegend'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot expected convergence rate
% noise-to-signal ratio
v = v_curve_array;

% solve the plot by taking derivative
% opt. step size
if(n==10)
    d =sprintf('n \\rightarrow \\infty');
    plot(v,c_mu_lambda*mu./(sqrt(1+v.*v)),'DisplayName',d);
    xlabel('noise-to-signal ratio \upsilon','FontSize',15);%
    ylabel('opt. normalized step size \sigma^*','FontSize',15); 
    set(gca, 'XScale', 'log');
    set(gca,'FontSize',15);
    p1 = sprintf('opt. normalized step size (%d/%d,%d)-ES',mu,mu,lambda);
    title(p1,'fontsize',20);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot convergence rate
figure(FIG_NUM);
hold on
d =sprintf('n = %.0f',n);
scatter(v_array, v_convergence_rate_array,typeDot,scatterColour,'DisplayName',d); hold on; 

legend('-DynamicLegend'); 
legend('show');
leg = legend();
leg.FontSize = 10;
title(leg,'dimension of data');
hold off
ylabel('convergnece rate c','fontsize',20);
xlabel('noise-to-signal ratio \upsilon','fontsize',20);
%set(gca,'xscale','log')
set(gca,'FontSize',15);
ylim([0,inf])                                                               % y starts from 0
d = sprintf("convergence rate (%d/%d,%d)-ES",mu,mu,lambda);
title(d,'FontSize', 25);


% save all fig
p2 = sprintf('opt_step_size_%d_%d_%d_ES.fig',mu,mu,lambda);
saveas(gcf,p2);

end



