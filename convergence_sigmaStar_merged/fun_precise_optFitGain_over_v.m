% add different curve for expected progress rate over dim n
% plot opt. normalized step size & opt. expected fitness gain over v 

function fun_precise_optFitGain_over_v(f,NUM_OF_RUNS,mu,lambda,v_array,n,scatterColour,typeDot,FIG_NUM,c_mu_lambda,v_curve_array)
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



% f = @(x) (x'*x);

%n = 10;                 % dim of the data
%ita = 4;     
% x0 = randn(n,mu); 

NUM_OF_ITERATIONS = 2000;
% NUM_OF_RUNS = 1;



% matrix for different sigma* matrix(sigma_i,median of NUM_OF_RUNS)
s_start = 0.2;
increment =0.2;
s_end = 32.2;

sigma_success_rate_array = zeros(1,(s_end-s_start)/increment+1);
sigma_T_array = zeros(1,(s_end-s_start)/increment+1);
sigma_f_x_array = zeros((s_end-s_start)/increment+1,6000);
sigma_counvergence_rate_array = zeros(1,(s_end-s_start)/increment+1);

% temp_sigma_success_rate_array = zeros(1,NUM_OF_RUNS);
% temp_sigma_T_array = zeros(1,NUM_OF_RUNS);
% temp_sigma_f_x_array = zeros(NUM_OF_RUNS,6000);
% temp_sigma_counvergence_rate_array = zeros(1,NUM_OF_RUNS);
% 


sigma_array = 0.00001*2.^(1:1:14);
SIGMA_LENGTH = length(sigma_array);

sigma_success_rate_array = zeros(1,SIGMA_LENGTH);
sigma_T_array = zeros(1,SIGMA_LENGTH);
sigma_f_x_array = zeros(SIGMA_LENGTH,6000);
sigma_counvergence_rate_array = zeros(1,SIGMA_LENGTH);

v_convergence_rate_array = zeros(1,V_LENGTH);                               % store the max median c over v
index_convergence_rate_array = zeros(1,V_LENGTH);                               % store the max median c over v

for k = 1:1:V_LENGTH 
    v_temp = v_array(k); 
    i = 1;
    for sigma_star = s_start:increment:s_end
%         k = 1:1:SIGMA_LENGTH
%         sigma_star = sigma_array(k);
        for j = 1:1:NUM_OF_RUNS
            x0 = randn(n,mu); 
            sigma_ep_star = v_temp*sigma_star;
            a = mml_noise(f,x0,sigma_star,sigma_ep_star,lambda,NUM_OF_ITERATIONS);
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
    [v_convergence_rate_array(k) index_convergence_rate_array(k)] = max(sigma_counvergence_rate_array);
    % STORE MAX index 
end

figure(FIG_NUM);
hold on
legend('-DynamicLegend'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot expected convergence rate for n-> infty
% noise-to-signal ratio
v = v_array;

% solve the plot by taking derivative
% opt. step size
if(n==100)
    d =sprintf('n \\rightarrow \\infty');
    plot(v,c_mu_lambda*mu./(sqrt(1+v.*v)),'Color','k','DisplayName',d);
    
    xlabel('noise-to-signal ratio \vartheta','FontSize',15);%
    ylabel('opt. normalized step size \sigma*','FontSize',15); 
    set(gca, 'XScale', 'log');
    set(gca,'FontSize',15);
    p1 = sprintf('opt. normalized step size (%d/%d,%d)-ES',mu,mu,lambda);
    title(p1,'fontsize',20);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot expected fitness gain for diffeent noise t
figure(FIG_NUM);
hold on;
sigma_star = 0.01:0.001:8;
% sigma_star = sigma_array;
v_trans = transpose(v); 
% a matrix 
% row: different noise-to-signal ratio v
% col: different normalized step size sigmaStar
expected_curve = c_mu_lambda*sigma_star.*(1+sigma_star.^2/2/mu/n)./(sqrt(1+sigma_star.^2/mu/n).*sqrt(1+v_trans.^2+sigma_star.^2/2/n))-n*(sqrt(1+sigma_star.^2/mu/n)-1);
[max_fitness_array s_index] = max(expected_curve,[],2);
max_sigma_array = zeros(length(v),1);
for i = 1:1:length(s_index)
    max_sigma_array(i) = sigma_star(s_index(i));
end
d1 = sprintf('n = %d',n);                                                % legend specify dim 
% estimation for specific n eq. (11) https://core.ac.uk/download/pdf/81976199.pdf
if(n==10) 
    plot(v,max_sigma_array,'-','Color',scatterColour,'DisplayName',d1);
elseif(n==100)
    plot(v,max_sigma_array,'-','Color',scatterColour,'DisplayName',d1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot convergence rate
figure(FIG_NUM);
hold on
d =sprintf('n = %.0f',n);
temp_sStar_array = s_start:increment:s_end;
sigma_max_exp = zeros(1,length(index_convergence_rate_array));
for i = 1:1:length(index_convergence_rate_array)
    sigma_max_exp(i) = temp_sStar_array(index_convergence_rate_array(i));
end
% scatter(v_array, v_convergence_rate_array,typeDot,scatterColour,'DisplayName',d); hold on; 
scatter(v_array, sigma_max_exp,typeDot,scatterColour,'DisplayName',d); hold on; 

legend('-DynamicLegend'); 
legend('show');
leg = legend();
leg.FontSize = 10;
title(leg,'dimension of data');
% ylabel('convergnece rate c','fontsize',20);
xlabel('noise-to-signal ratio \vartheta','fontsize',20);
ylabel('opt. normalized step size \sigma*','FontSize',15); 
%set(gca,'xscale','log')
set(gca,'FontSize',15);
ylim([0,inf])                                                               % y starts from 0
d = sprintf("convergence rate (%d/%d,%d)-ES",mu,mu,lambda);
title(d,'FontSize', 20);


% save all fig
p2 = sprintf('1opt_normalized_step_size_%d_%d_%d_ES.fig',mu,mu,lambda);
saveas(gcf,p2);


end



