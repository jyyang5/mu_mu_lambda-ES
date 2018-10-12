% add expected c (curve), mml (dotted line) 
% plot convergence rate and success rate 
% for different sigma* and v = sigma_ep*/sigma*

function fun_fitness_sigmaStar_multi(NUM_OF_RUNS,mu,lambda,ita,n,scatterColour,typeDot,LINE_OR_NOT,c_mu_lambda)
%Input
%   ita:                  % ita = sigma_ep_star /sigma_star
%   n:                    dim of data
%   scatterColour:        colour of the plots
%   typeDot:              type of dot (+  Plus sign, o  Circle, *  Asterisk, x  Cross, . dotted line)
%   LINE_OR_NOT:          justScatter = 0, includeLine = 1
%   c_mu_lambda:          a factor for computing n-> infty convergence rate
%Return
%   scatter or scatter+curve

f = @(x) (x'*x);

%n = 10;                 % dim of the data
%ita = 4;     
x0 = randn(n,mu); 

NUM_OF_ITERATIONS = 2000;
% NUM_OF_RUNS = 2;




% matrix for different sigma* matrix(sigma_i,median of NUM_OF_RUNS)
s_start = 0.2;
increment =0.2;
s_end = 8.2;
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
        sigma_ep_star = ita*sigma_star;
        if(typeDot == '.')
            a = mml_normal(f,x0,sigma_star,sigma_ep_star,lambda,NUM_OF_ITERATIONS);
        else 
            a = mml_noise(f,x0,sigma_star,sigma_ep_star,lambda,NUM_OF_ITERATIONS);
        end
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
    %d = sprintf("sigma_star = %d -> success rate = %d", sigma_star,cell2mat(a(5)));
    %disp(d);
    
end

% % success rate
% %figure(1)
% hold on
% scatter(s_start:increment:(s_end), sigma_success_rate_array,typeDot,scatterColour);
% 
% % if n -> \infty plot the curve
% if LINE_OR_NOT==1
%     plot(s_start:increment:(s_end), sigma_success_rate_array,scatterColour);
% end
% % each ita plot ita = some n
% if n ==10
%     d = sprintf('v = %d', ita);
%     text((s_start+s_end)/2, median(sigma_success_rate_array), d, 'Color', scatterColour);
% end
% 
% hold off
% ylabel('sucess rate c');
% xlabel('normalized step size \sigma^*');
% %set(gca,'xscale','log')
% d = sprintf("success rate over log(sigma*) with sigma_{ep} = %d",sigma_ep_star);
% title(d);
% 
% end

% convergence rate
%figure(1)
figure(4);
legend('-DynamicLegend'); 
hold on
% dotted line (normal mml-ES)
if(typeDot == '.')
    d = sprintf('without GP');
    plot(s_start:increment:(s_end), sigma_counvergence_rate_array,':','DisplayName',d); hold on; 
% need scatter plots (through real run)
else
    d = sprintf('v = %.2f',ita);
    scatter(s_start:increment:(s_end), sigma_counvergence_rate_array,typeDot,scatterColour,'DisplayName',d); hold on; 
end
% if n -> \infty plot the curve
if LINE_OR_NOT==1
    sigma_star = s_start:increment:s_end;
    expected_cure = sigma_star*(c_mu_lambda)./sqrt(1+ita^2)-sigma_star.^2./(2*mu);
    plot(sigma_star,expected_cure,'Color',scatterColour,'DisplayName',d); hold on;    
%     plot(s_start:increment:(s_end), sigma_counvergence_rate_array,scatterColour);
end
% each ita plot ita = some n
if n ==10
%     d = sprintf('\nu = %d', ita);
    text((s_start+s_end)/2-1, median(sigma_counvergence_rate_array), d, 'Color', scatterColour, 'fontsize',15);
end
legend('-DynamicLegend'); 
legend('show');
leg = legend();
leg.FontSize = 10;
title(leg,'noise-to-signal-ratio \upsilon');
hold off
ylabel('convergnece rate c','fontsize',20);
xlabel('normalized step size \sigma^*','fontsize',20);
%set(gca,'xscale','log')
set(gca,'FontSize',15);
ylim([0,inf])
d = sprintf("convergence rate (%d/%d,%d)-ES",mu,mu,lambda);
title(d,'FontSize', 25);
p2 = sprintf('fitGain_%d_%d_%d_ES.fig',mu,mu,lambda);
saveas(gcf,p2);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % convergence rate
% subplot(1,3,2);
% plot(s_start:increment:(s_end), sigma_counvergence_rate_array);
% ylabel('convergence rate');
% xlabel('log(sigma*)');
% set(gca,'xscale','log')
% d = sprintf("convergence rate over log(sigma*) with sigma_{ep} = %d",sigma_ep_star);
% title(d);
% % number of objective function evaluations
% subplot(1,3,3);
% plot(s_start:increment:(s_end), sigma_T_array);
% ylabel('number of objective function evaluation');
% xlabel('log(sigma*)');
% set(gca,'xscale','log')
% d = sprintf("number of objective function evaluation over log(sigma*) with sigma_{ep} = %d",sigma_ep);
% title(d);
% 
% d0 = sprintf('Proportion of parent solution not the lowest for the quadratic model = %d',median(parent_not_lowest_of_quadratic_array));
% disp(d0);


% 
% % matrix for different c matrix(c_i,median of NUM_OF_RUNS)
% sep_start = -10; 
% sep_end = 3;
% sep_num = sep_end-sep_start+1;
% sep_range = 10.^[sep_start:sep_end];
% sep_success_rate_matrix = zeros(NUM_OF_RUNS,sep_num);
% sep_counvergence_rate_matrix = zeros(NUM_OF_RUNS,sep_num);
% sep_T_array = zeros(NUM_OF_RUNS,sep_num);
% 
% 
% % over c
% temp_c_success_rate_array = zeros(1,NUM_OF_RUNS);
% temp_sep_T_array = zeros(1,NUM_OF_RUNS);
% temp_c_f_x_array = zeros(NUM_OF_RUNS,6000);
% temp_c_counvergence_rate_array = zeros(1,NUM_OF_RUNS);
% 
% sep_success_rate_array = zeros(1,sep_num);
% sep_T_array = zeros(1,sep_num);
% sep_f_x_array = zeros(sep_num,6000);
% sep_counvergence_rate_array = zeros(1,sep_num);
% 
% for i = 1:1:sep_num
%     sigma_exp_temp = sep_range(i);
%     for j = 1:1:NUM_OF_RUNS
%         a = strategy(f,x0,sigma0,sigma_exp_temp,NUM_OF_ITERATIONS);
%         temp_sep_f_x_array(j,:) = cell2mat(a(3));
%         % avoid nan in fx (if nan appears -> set 0)
%         if(sum(isnan(temp_sep_f_x_array(j,:))) == 0)
%             temp_sep_counvergence_rate_array(j) = cell2mat(a(7));
%             temp_sep_T_array(j) = cell2mat(a(1));
%             temp_sep_success_rate_array(j) = cell2mat(a(5));
%         else 
%             temp_sep_counvergence_rate_array(j) = 0;
%             temp_sep_T_array(j) = 0;
%             temp_sep_success_rate_array(j) = 0;
%         end
%     end
%     for j = 1:1:6000
%         sep_f_x_array(i,j) = median(temp_sep_f_x_array(:,j));
%     end
%     sep_success_rate_array(i) = median(temp_sep_success_rate_array);
%     sep_counvergence_rate_array(i) = median(temp_sep_counvergence_rate_array);
%     sep_T_array(i) = median(temp_sep_T_array);
%     
%     
%     i = i + 1;
%     %d = sprintf("sigma_star = %d -> success rate = %d", sigma_star,cell2mat(a(5)));
%     %disp(d);
%     
% end
% 
% % success rate
% figure(3)
% subplot(3,1,1);
% plot(sep_range, sep_success_rate_array);
% ylabel('sucess rate');
% xlabel('log(sigma_ep)');
% set(gca,'xscale','log')
% d = sprintf('success rate over noise rate log(sigma_{ep}) with sigma* = %d',sigma0);
% title('success rate over noise rate log(sigma_{ep})');
% % convergence rate
% subplot(3,1,2);
% plot(sep_range, sep_counvergence_rate_array);
% ylabel('convergence rate');
% xlabel('log(sigma_{ep})');
% set(gca,'xscale','log')
% d = sprintf("convergence rate over noise rate log(sigma_{ep}) with sigma* = %d",sigma0)
% title(d);
% % number of objective function evaluations
% subplot(3,1,3);
% plot(sep_range, sep_T_array);
% ylabel('number of objective function evaluation');
% xlabel('log(sigma_{ep})');
% set(gca,'xscale','log')
% d = sprintf("number of objective function evaluation over noise rate log(sigma_{ep}) with sigma* = %d",sigma0)
% title(d);
% 
