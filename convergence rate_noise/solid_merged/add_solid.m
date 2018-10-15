% Plot solid line in Expected fitness gain over sigma* and sigma_ep^* 
% solid line in convergence rate graph

% a = openfig('test.fig');

mu = 3;
lambda = 10;                            % not used but for computing c_mu_lambda

step = 0.0000001;
x = -10:step:10;
% expected fitness gain for n->\infty  under mu and lambda
% c_mu_lambda = (lambda-mu)/(2*pi)*nchoosek(lambda,mu)*sum(exp(-x.^2).*(normcdf(x)).^(lambda-mu-1).*(1-normcdf(x)).^(mu-1))*step;
c_mu_lambda = 1.065389626877247; % expected convergence rate for (3/3,10)-ES

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
figure(1);
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opt. expected progress rate and normalized step size calculated
% numerically 

% noise-to-signal ratio
s_start = 0.01;
increment =0.001;
s_end = 8;
sigma_star = s_start:increment:s_end;
v = transpose(0.1:0.001:50); 
% a matrix 
% row: different noise-to-signal ratio v
% col: different normalized step size sigmaStar
n =10;
expected_cure = c_mu_lambda*sigma_star.*(1+sigma_star.^2/2/mu/n)./(sqrt(1+sigma_star.^2/mu/n).*sqrt(1+v.^2+sigma_star.^2/2/n))-n*(sqrt(1+sigma_star.^2/mu/n)-1);
[max_fitness_array s_index] = max(expected_cure,[],2);
max_sigma_array = zeros(length(v),1);
for i = 1:1:length(s_index)
    max_sigma_array(i) = sigma_star(s_index(i));
end

n =100;
expected_cure1 = c_mu_lambda*sigma_star.*(1+sigma_star.^2/2/mu/n)./(sqrt(1+sigma_star.^2/mu/n).*sqrt(1+v.^2+sigma_star.^2/2/n))-n*(sqrt(1+sigma_star.^2/mu/n)-1);
[max_fitness_array1 s_index1] = max(expected_cure1,[],2);
max_sigma_array1 = zeros(length(v),1);
for i = 1:1:length(s_index1)
    max_sigma_array1(i) = sigma_star(s_index1(i));
end



v = transpose(v);
% v = 0.1:0.00001:10;

% solve the plot by taking derivative
% opt. step size
figure(2);
hold on;
% precise
plot(v,max_sigma_array,'r');
plot(v,max_sigma_array1,'b');
% original
plot(v,c_mu_lambda*mu./(sqrt(1+v.*v)),'k');
hold off;
xlabel('noise-to-signal ratio \upsilon','FontSize',15);%
ylabel('opt. normalized step size \sigma^*','FontSize',15); 
legend('n=10','n=100','n \rightarrow \infty');
set(gca, 'XScale', 'log');
set(gca,'FontSize',15);



% opt. fitness gain given opt. step size
figure(3);
% precise
hold on;
plot(v,max_fitness_array,'r');
plot(v,max_fitness_array1,'b');
% original
sigma_star = c_mu_lambda*mu./(sqrt(1+v.*v));
plot(v,sigma_star*(c_mu_lambda)./sqrt(1+v.^2)-sigma_star.^2./(2*mu),'k');
hold off;
ylabel('opt. expected fitness gain \eta','FontSize',15);
xlabel('noise-to-signal-ratio \upsilon','FontSize',15); 
legend('n=10','n=100','n \rightarrow \infty');
set(gca, 'XScale', 'log');
set(gca,'FontSize',15);





