% plot mml-ES with GP over sigmaStar
%

function fun_GP_fitness_sigmaStar(f,NUM_OF_RUNS,mu,lambda,n,scatterColour,FIGURE_NUM)
%Input
%   n:                    dim of data (n==10 dashed line, n==100 dotted line)
%   scatterColour:        colour of the plots
%   LINE_OR_NOT:          justScatter = 0, includeLine = 1
%   c_mu_lambda:          a factor for computing n-> infty convergence rate
%   FIGURE_NUM:           # figure
%Return
%   plot using GP

% f = @(x) (x'*x);

%n = 10;                 % dim of the data
%ita = 4;     
x0 = randn(n,mu); 

NUM_OF_ITERATIONS = 2000;
TRAINING_SIZE = 40;
% NUM_OF_RUNS = 2;




% matrix for different sigma* matrix(sigma_i,median of NUM_OF_RUNS)
s_start = 0.2;
% increment =0.5;
if(lambda<=10)
    s_end = 6.2;
elseif(lambda<=20)
    s_end = 8.2;
else 
    s_end = 10.2;
end
increment =0.2;
s_end = 7.2;
% s_start = 1;
% increment =1;
% s_end = 8;
temp_sigma_success_rate_array = zeros(1,NUM_OF_RUNS);
temp_parent_not_lowest_of_quadratic_array = zeros(1,NUM_OF_RUNS);
temp_sigma_T_array = zeros(1,NUM_OF_RUNS);
temp_sigma_f_x_array = zeros(NUM_OF_RUNS,6000);
temp_sigma_counvergence_rate_array = zeros(1,NUM_OF_RUNS);

sigma_success_rate_array = zeros(1,uint8((s_end-s_start)/increment+1));
parent_not_lowest_of_quadratic_array = zeros(1,uint8((s_end-s_start)/increment+1));
sigma_T_array = zeros(1,uint8((s_end-s_start)/increment+1));
sigma_f_x_array = zeros(uint8((s_end-s_start)/increment+1),6000);
sigma_counvergence_rate_array = zeros(1,uint8((s_end-s_start)/increment+1));

i = 1;
for sigma_star = s_start:increment:s_end
    for j = 1:1:NUM_OF_RUNS
        x0 = randn(n,mu); 
        % GP mml-ES
        a = mml_GP(f,x0,sigma_star,TRAINING_SIZE,lambda,NUM_OF_ITERATIONS);    
        %temp_parent_not_lowest_of_quadratic_array(j) = cell2mat(a(6));
        %temp_sigma_counvergence_rate_array(j) = cell2mat(a(7));
        %temp_sigma_T_array(j) = cell2mat(a(1));
        %temp_sigma_success_rate_array(j) = cell2mat(a(5));
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
    for j = 1:1:6000
        sigma_f_x_array(i,j) = median(temp_sigma_f_x_array(:,j));
    end
    sigma_success_rate_array(i) = median(temp_sigma_success_rate_array);
    parent_not_lowest_of_quadratic_array(i) = median(temp_parent_not_lowest_of_quadratic_array);
    sigma_counvergence_rate_array(i) = median(temp_sigma_counvergence_rate_array);
    sigma_T_array(i) = median(temp_sigma_T_array);
    
    
    i = i + 1;
    
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot progress rate for GP mml-ES
figure(FIGURE_NUM);
legend('-DynamicLegend'); 
hold on
d = sprintf('withGP n=%d',n);

if(n==10)
    plot(s_start:increment:s_end, sigma_counvergence_rate_array,'--','Color',scatterColour,'DisplayName',d); hold on; 
elseif(n==100)
    plot(s_start:increment:s_end, sigma_counvergence_rate_array,':','Color',scatterColour,'DisplayName',d); hold on; 
end
% each ita plot theta = some n
% if n ==10
% %     d = sprintf('\nu = %d', ita);
%     text((s_start+s_end)/2-1, median(sigma_counvergence_rate_array), d, 'Color', scatterColour, 'fontsize',15);
% end

legend('-DynamicLegend'); 
legend('show');
leg = legend();
leg.FontSize = 10;
title(leg,'noise-to-signal-ratio \vartheta');
hold off
ylabel('expected fitness gain \eta','fontsize',20);
xlabel('normalized step size \sigma^*','fontsize',20);
%set(gca,'xscale','log')
set(gca,'FontSize',15);
ylim([0,inf])
d = sprintf("expected normalized step size (%d/%d,%d)-ES",mu,mu,lambda);
title(d,'FontSize', 20);
p2 = sprintf('1fitGain_%d_%d_%d_ES.fig',mu,mu,lambda);
saveas(gcf,p2);

end

