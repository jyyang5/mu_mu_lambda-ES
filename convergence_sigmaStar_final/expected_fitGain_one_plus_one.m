% Plot dotted line in Expected fitness gain over sigma* for [(1+1)-ES]
% solid line in convergence rate graph

% a = openfig('test.fig');
%(3/3,10),(5/5,20), (10,40)
mu = 10;
lambda = 40;                            % not used but for computing c_mu_lambda


% First attempt 
% sigma_star_temp = transpose(1:0.1:5);
% z_end = 10;
% z_step = 0.00001;
% % only one sigma
% z = sigma_star_temp/2:z_step:z_end;
% % Probablity of making an objective function call
% p_eval = normcdf(-sigma_star.^2/.2./(sqrt(sigma_star.^2+sigma_ep_star.^2)));
% 
% % Expected value of change in objective function 
% sigma_star = sigma_star_temp;
% fz = (sigma_star*z-sigma_star.^2/2)*exp(-z.^2/2)*normcdf((sigma_star.*z-sigma_star.^2./2)./sigma_ep_star);
% expected_delta = 1/sqrt(2*pi)*sum(fz)*z_step;
% 
% vartheta = expected_delta./p_eval;

% Iterate over v

v = 0.1;
sigma_star_array = transpose(0.1:0.01:1.8);
z_start = 0.001;
z_step = 0.001; 
z_end = 10;
z_array = z_start:z_step:z_end;
z_LENGTH = length(z_array);
z = z_array+sigma_star_array;

% sigma_star_temp = transpose(1:0.1:5);
% sigma_star_start = 0.2;
% z_step = 0.01;
% z_end = 10;

% only one sigma
% z = sigma_star_temp/2:z_step:z_end;
% Probablity of making an objective function call
p_eval = normcdf(-1./2./sqrt(v.^2));

% Expected value of change in objective function 
sigma_star = sigma_star_array;
fz = (sigma_star*z-sigma_star.^2/2)*exp(-z.^2/2)*normcdf(v.*(z-transpose(repmat(sigma_star,1,z_LENGTH)/2)));
expected_delta = 1/sqrt(2*pi)*sum(fz)*z_step;

vartheta = expected_delta./p_eval;



% 
% % expected fitness gain for n->\infty  under mu and lambda
% c_mu_lambda = (lambda-mu)/(2*pi)*nchoosek(lambda,mu)*sum(exp(-x.^2).*(normcdf(x)).^(lambda-mu-1).*(1-normcdf(x)).^(mu-1))*step;
% % c_mu_lambda = 1.065389626877247; % expected convergence rate for (3/3,10)-ES
% % c_mu_lambda = 1.214478382788638; % expected convergence rate for (5/5,20)-ES
% % c_mu_lambda = 1.242204493664515; % expected convergence rate for (10/10,40)-ES
% 
% 
% % range for sigma_satr
% s_start = 0.2;
% increment =0.2;
% s_end = 6.2;
% sigma_star = s_start:increment:s_end;
% 
% % noise-to-signal ratio = 4
% ita = 4;
% final = sigma_star*(c_mu_lambda)./sqrt(1+ita^2)-sigma_star.^2./(2*mu);
% 
% % noise-to-signal ratio = 1
% ita = 1;
% final1 = sigma_star*(c_mu_lambda)./sqrt(1+ita^2)-sigma_star.^2./(2*mu);
% % noise-to-signal ratio = 1/4
% ita = 1/4;
% final2 = sigma_star*(c_mu_lambda)./sqrt(1+ita^2)-sigma_star.^2./(2*mu);
% % noise-to-signal ratio = 0
% ita = 0;
% final3 = sigma_star*(c_mu_lambda)./sqrt(1+ita^2)-sigma_star.^2./(2*mu);
% 
% % expected fitness gain over noise-to-signal ratio v
% figure(2);
% hold on;
% plot(sigma_star,final,'g');     % ita = 4
% plot(sigma_star,final1,'r');    % ita = 1
% plot(sigma_star,final2,'b');    % ita = 1/4
% plot(sigma_star,final3,'k');    % ita = 0
% hold off;
% xlabel('\sigma*','FontSize',15);%
% ylabel('Expected fitness gain','FontSize',15); 
% ylim([0 inf]);
% set(gca,'FontSize',15);
% 
% 
% 
% % noise-to-signal ratio
% v = 0.1:0.00001:10;
% 
% % solve the plot by taking derivative
% % opt. step size
% figure(3);
% plot(v,c_mu_lambda*mu./(sqrt(1+v.*v)));
% xlabel('noise-to-signal ratio \upsilon','FontSize',15);%
% ylabel('opt. normalized step size \sigma^*','FontSize',15); 
% set(gca, 'XScale', 'log');
% set(gca,'FontSize',15);
% 
% 
% 
% % opt. fitness gain given opt. step size
% figure(4);
% sigma_star = c_mu_lambda*mu./(sqrt(1+v.*v));
% plot(v,sigma_star*(c_mu_lambda)./sqrt(1+v.^2)-sigma_star.^2./(2*mu));
% ylabel('opt. expected fitness gain \eta','FontSize',15);
% xlabel('noise-to-signal-ratio \upsilon','FontSize',15); 
% set(gca, 'XScale', 'log');
% set(gca,'FontSize',15);
% 
% 
% 
% 
% 