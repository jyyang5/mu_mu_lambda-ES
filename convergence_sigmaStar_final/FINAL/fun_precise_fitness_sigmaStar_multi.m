% use new n-> infty for expected convergence rate over different dim 
% add expected fitness gain n-> infty(curve), mml (dotted line), scatter n= 10,100 
% plot convergence rate and success rate 
% for different sigma* and v = sigma_ep*/sigma*

function fun_precise_fitness_sigmaStar_multi(f,NUM_OF_RUNS,mu,lambda,ita,n,scatterColour,typeDot,LINE_OR_NOT,c_mu_lambda) 
%Input
%   ita:                  noise-to-signal ratio: ita = sigma_ep_star /sigma_star
%   n:                    dim of data
%   scatterColour:        colour of the plots
%   typeDot:              type of dot (+  Plus sign, o  Circle, *  Asterisk, x  Cross, . dotted line)
%                         normal mml-ES = '.'
%                         GP mml-ES = '-'  
%   LINE_OR_NOT:          justScatter = 0, includeLine = 1
%   c_mu_lambda:          a factor for computing n-> infty convergence rate
%Return
%   scatter or scatter+curve

% f = @(x) (x'*x);

%n = 10;                 % dim of the data
%ita = 4;     
x0 = randn(n,mu); 

NUM_OF_ITERATIONS = 2000;
% NUM_OF_RUNS = 2;




% matrix for different sigma* matrix(sigma_i,median of NUM_OF_RUNS)
s_start = 0.0001;
increment =0.2;
s_end = 10+s_start;
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
        sigma_ep_star = ita*sigma_star;
        % normal mml-ES
        if(typeDot == '.')
            a = mml_normal(f,x0,sigma_star,sigma_ep_star,lambda,NUM_OF_ITERATIONS);
        else 
            a = mml_noise(f,x0,sigma_star,sigma_ep_star,lambda,NUM_OF_ITERATIONS);
        end
        
        temp_sigma_f_x_array(j,:) = cell2mat(a(3));
        % avoid nan in fx (if nan appears -> set 0)
        if(sum(isnan(temp_sigma_f_x_array(j,:))) == 0)
            temp_sigma_counvergence_rate_array(j) = cell2mat(a(7));
            temp_sigma_T_array(j) = cell2mat(a(1));
        else 
            temp_sigma_counvergence_rate_array(j) = 0;
            temp_sigma_T_array(j) = 0;
        end
    end
    for j = 1:1:6000
        sigma_f_x_array(i,j) = median(temp_sigma_f_x_array(:,j));
    end
    sigma_success_rate_array(i) = median(temp_sigma_success_rate_array);
    parent_not_lowest_of_quadratic_array(i) = median(temp_parent_not_lowest_of_quadratic_array);
    sigma_counvergence_rate_array(i) = median(temp_sigma_counvergence_rate_array)./lambda;
    sigma_T_array(i) = median(temp_sigma_T_array);    
    i = i + 1;

    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot precise expected progress rate when n -> infty  
figure(4);
legend('-DynamicLegend'); 
hold on
% need scatter plots (through experiemnt runs different noise-to signal ratio and dim)
scatter_range=(1:5:(s_end-s_start)/increment+1);
sigma_temp_array = s_start:increment:s_end;
d = sprintf('\\vartheta = %.2f N=%d',ita,n);
scatter(sigma_temp_array(scatter_range), sigma_counvergence_rate_array(scatter_range),typeDot,scatterColour,'DisplayName',d); hold on; 

% if n -> \infty plot the curve(expected fitness gain)
if LINE_OR_NOT==1
    hold on;
    % first approximation n->infty
    sigma_star = s_start:increment:s_end;
    expected_cure = sigma_star*(c_mu_lambda)./sqrt(1+ita^2)-sigma_star.^2./(2*mu);
    d1 = sprintf('\\vartheta = %.2f N \\rightarrow \\infty',ita);
    plot(sigma_star,expected_cure./lambda,'Color',scatterColour,'DisplayName',d1); hold on; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot curve for mml-ES without model assistance [simply /lambda]
    if(ita==0)   
        d5 = sprintf('No model N \\rightarrow \\infty');
        plot(sigma_star,expected_cure./lambda./lambda,':','Color','k','DisplayName',d5); hold on;
        d5 = sprintf('No model');
        text(5, max(sigma_counvergence_rate_array./lambda),d5, 'Color', 'k', 'fontsize',15);
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot precise expected progress rate with a specified dimension n 
% d1 = sprintf('\\vartheta = %.2f n=%d',ita,n);
% v = ita;
% hold on;
% sigma_star = s_start:increment:s_end;
% precise_curve = c_mu_lambda*sigma_star.*(1+sigma_star.^2/2/mu/n)./(sqrt(1+sigma_star.^2/mu/n).*sqrt(1+v.^2+sigma_star.^2/2/n))-n*(sqrt(1+sigma_star.^2/mu/n)-1);
% if(n==10 && typeDot~='.')
%     plot(sigma_star,precise_curve,'--','Color',scatterColour,'DisplayName',d1);
% elseif(n==100 && typeDot~='.')
%     plot(sigma_star,precise_curve,':','Color',scatterColour,'DisplayName',d1);
% end

% each ita plot theta = some n
if n ==10
    d = sprintf('\\vartheta = %.2f', ita);
    text(5, max(sigma_counvergence_rate_array), d, 'Color', scatterColour, 'fontsize',15);
end


legend('-DynamicLegend'); 
legend('show');
leg = legend();
leg.FontSize = 10;
% title(leg,'Noise-to-signal-ratio  Dimension');
hold off
ylabel('expected fitness gain \eta','fontsize',20);
xlabel('normalized step size \sigma^*','fontsize',20);
%set(gca,'xscale','log')
set(gca,'FontSize',15);
ylim([0,inf])
xlim([0 10]);
box on;

% d = sprintf("expected normalized step size (%d/%d,%d)-ES",mu,mu,lambda);
% title(d,'FontSize', 20);
% p2 = sprintf('1fitGain_%d_%d_%d_ES.fig',mu,mu,lambda);
% saveas(gcf,p2);

end

