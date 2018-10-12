% Plot solid line in Expected fitness gain over sigma* and sigma_ep^* 
% solid line in convergence rate graph

% a = openfig('test.fig');
%(3/3,10),(5/5,20), (10,40)
mu = 10;
lambda = 40;                            % not used but for computing c_mu_lambda

step = 0.0000001;
x = -10:step:10;
% expected fitness gain for n->\infty  under mu and lambda
c_mu_lambda = (lambda-mu)/(2*pi)*nchoosek(lambda,mu)*sum(exp(-x.^2).*(normcdf(x)).^(lambda-mu-1).*(1-normcdf(x)).^(mu-1))*step;
% c_mu_lambda = 1.065389626877247; % expected convergence rate for (3/3,10)-ES
% c_mu_lambda = 1.214478382788638; % expected convergence rate for (5/5,20)-ES
% c_mu_lambda = 1.242204493664515; % expected convergence rate for (10/10,40)-ES


% range for sigma_satr
s_start = 0.2;
increment =0.2;
s_end = 6.2;
sigma_star = s_start:increment:s_end;

% noise-to-signal ratio = 4
ita = 4;
final = sigma_star*(c_mu_lambda)./sqrt(1+ita^2)-sigma_star.^2./(2*mu);

% noise-to-signal ratio = 1
ita = 1;
final1 = sigma_star*(c_mu_lambda)./sqrt(1+ita^2)-sigma_star.^2./(2*mu);
% noise-to-signal ratio = 1/4
ita = 1/4;
final2 = sigma_star*(c_mu_lambda)./sqrt(1+ita^2)-sigma_star.^2./(2*mu);
% noise-to-signal ratio = 0
ita = 0;
final3 = sigma_star*(c_mu_lambda)./sqrt(1+ita^2)-sigma_star.^2./(2*mu);

% expected fitness gain over noise-to-signal ratio v
figure(2);
hold on;
plot(sigma_star,final,'g');     % ita = 4
plot(sigma_star,final1,'r');    % ita = 1
plot(sigma_star,final2,'b');    % ita = 1/4
plot(sigma_star,final3,'k');    % ita = 0
hold off;
xlabel('\sigma*','FontSize',15);%
ylabel('Expected fitness gain','FontSize',15); 
ylim([0 inf]);
set(gca,'FontSize',15);



% noise-to-signal ratio
v = 0.1:0.00001:10;

% solve the plot by taking derivative
% opt. step size
figure(3);
plot(v,c_mu_lambda*mu./(sqrt(1+v.*v)));
xlabel('noise-to-signal ratio \upsilon','FontSize',15);%
ylabel('opt. normalized step size \sigma^*','FontSize',15); 
set(gca, 'XScale', 'log');
set(gca,'FontSize',15);



% opt. fitness gain given opt. step size
figure(4);
sigma_star = c_mu_lambda*mu./(sqrt(1+v.*v));
plot(v,sigma_star*(c_mu_lambda)./sqrt(1+v.^2)-sigma_star.^2./(2*mu));
ylabel('opt. expected fitness gain \eta','FontSize',15);
xlabel('noise-to-signal-ratio \upsilon','FontSize',15); 
set(gca, 'XScale', 'log');
set(gca,'FontSize',15);





