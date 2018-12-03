% % Plot dotted line in Expected fitness gain over sigma* for [(1+1)-ES]
% % solid line in convergence rate graph
% 
% 
% mu = 10;
% lambda = 40;                            % not used but for computing c_mu_lambda
% 
% 
% v = 0.1;
% sigma_star_array = transpose(0.1:0.01:1.8);
% z_start = 0.001;
% z_step = 0.001; 
% z_end = 10;
% z_array = z_start:z_step:z_end;
% z_LENGTH = length(z_array);
% z = z_array+sigma_star_array;
% 
% 
% % Probablity of making an objective function call
% p_eval = normcdf(-1./2./sqrt(v.^2));
% 
% % Expected value of change in objective function 
% sigma_star = sigma_star_array;
% fz = (sigma_star*z-sigma_star.^2/2)*exp(-z.^2/2)*normcdf(v.*(z-transpose(repmat(sigma_star,1,z_LENGTH)/2)));
% expected_delta = 1/sqrt(2*pi)*sum(fz)*z_step;
% 
% vartheta = expected_delta./p_eval;
% 

% sigma_star_array = transpose(0.1:0.01:1.8);
% z_start = 0.001;
% z_step = 0.001; 
% z_end = 10;
% z_array = z_start:z_step:z_end;
% z_LENGTH = length(z_array);
% z = z_array+sigma_star_array;

sigma_star = 0.01:0.01:8.01;
sigma_star_trans = transpose(sigma_star);
z_start = 0;
z_step = 0.1; 
z_end = 20;
% z_step = (z_end-z_start)/z_LENGTH; 
z = (z_start:z_step:z_end)+sigma_star_trans./2;% sigma_star/2:(20-sigma_star/2)/z_LENGTH:20;
z_LENGTH = length(z_start:z_step:z_end);


% for v=0.1:0.1:10
% % Probablity of making an objective function call
% p_eval = normcdf(-sigma_star_trans./2./sqrt(1+v.^2));
% 
% fz = (repmat(sigma_star_trans,1,z_LENGTH).*z-repmat(sigma_star_trans.^2/2,1,z_LENGTH)).*exp(-z.^2/2).*normcdf(1./v.*(z-repmat(sigma_star_trans./2,1,z_LENGTH)));
% expected_delta = 1/sqrt(2*pi)*sum(fz,2)*z_step;
% result = expected_delta./p_eval;
% fprintf("v = %d, max=%d\n",v,max(result));
% end

v = 4;
p_eval = normcdf(-sigma_star_trans./2./sqrt(1+v.^2));
figure(1);
plot(sigma_star_trans,p_eval);

fz = (repmat(sigma_star_trans,1,z_LENGTH).*z-repmat(sigma_star_trans,1,z_LENGTH).^2/2).*exp(-z.^2/2).*normcdf(1./v.*(z-repmat(sigma_star_trans./2,1,z_LENGTH)));
expected_delta = 1/sqrt(2*pi)*sum(fz,2).*z_step;
result = expected_delta./p_eval;
figure(2);
plot(sigma_star_trans,result);





% fz = (sigma_star*z-sigma_star.^2/2)*exp(-z.^2/2)*normcdf(1./v.*(z_array-sigma_star.^2./2));
% expected_delta = 1/sqrt(2*pi)*sum(fz)*z_step;

% z_step = 0.0001;
v=1;
% Matrix approach
eta1 = @(sigma_star) 1/sqrt(2*pi)*sum((repmat(transpose(sigma_star),1,length(0:z_step:20)).*((0:z_step:20)+transpose(sigma_star./2))...
    -repmat(transpose(sigma_star),1,length(0:z_step:20)).^2/2).*exp(-((0:z_step:20)+...
    transpose(sigma_star)).^2/2).*normcdf(1./v.*(((0:z_step:20)+transpose(sigma_star./2))-...
    repmat(transpose(sigma_star),1,length(0:z_step:20))./2)),2).*z_step./normcdf(-transpose(sigma_star)./2./sqrt(1+v.^2));
% sigma_star = 1:2:5;
% sigma_star = 2;
v=1;
fminsearch(eta1,3)
% eta1 = @(sigma_star,v) 1/sqrt(2*pi)*sum((sigma_star.*z_array-sigma_star.^2/2)*exp(-z_array.^2/2)*normcdf(1./v.*(z_array-sigma_star.^2./2)))*z_step./normcdf(-1./2./sqrt(v.^2));
% temp = 1/sqrt(2*pi)*sum((sigma_star.*(sigma_star/2:(20-sigma_star/2)/z_LENGTH:20)-sigma_star.^2/2)*exp(-(sigma_star/2:(20-sigma_star/2)/z_LENGTH:20).^2/2)*normcdf(1./v.*((sigma_star/2:(20-sigma_star/2)/z_LENGTH:20)-sigma_star.^2./2)))*((20-sigma_star/2)/z_LENGTH)./normcdf(-1./2./sqrt(v.^2))
% eta1 = @(sigma_star,v) 1/sqrt(2*pi)*sum((sigma_star.*(sigma_star/2:(20-sigma_star/2)/z_LENGTH:20)-sigma_star.^2/2)*exp(-(sigma_star/2:(20-sigma_star/2)/z_LENGTH:20).^2/2)*normcdf(1./v.*((sigma_star/2:(20-sigma_star/2)/z_LENGTH:20)-sigma_star.^2./2)))*((20-sigma_star/2)/z_LENGTH)./normcdf(-1./2./sqrt(v.^2));


% sigma_star = [0,2];
% fminsearch(eta1,sigma_star)


function val = eta(sigma_star,v)
    z_start = sigma_star/2;
    z_LENGTH = 10000;
    z_end = 20;
    z_step = (z_end-z_start)/z_LENGTH; 
    z_array = z_start:z_step:z_end;
    % Probablity of making an objective function call
    p_eval = normcdf(-1./2./sqrt(v.^2));
    
    fz = (sigma_star*z-sigma_star.^2/2)*exp(-z.^2/2)*normcdf(1./v.*(z_array-sigma_star.^2./2));
    expected_delta = 1/sqrt(2*pi)*sum(fz)*z_step;
    val = expected_delta./p_eval;
    
end



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
