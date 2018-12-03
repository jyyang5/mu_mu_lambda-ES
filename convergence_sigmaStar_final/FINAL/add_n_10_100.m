
close all;
% Plot the experimental result (crosses and circles)
exp_v_array = v_array(1:10:101);
sigma_star_trans = opt_sigma_star_array(1:10:101);
sigma_star_array = transpose(sigma_star_trans);



% z = (z_start:z_step:z_end)+sigma_star_trans./2;                             % Put into a matrix strating from opt. sigma*/2
% z_LENGTH = length(z);
% p_eval = normcdf(-sigma_star_array./2./sqrt(1+exp_v_array.^2));
% v = exp_v_array;
% fz = (repmat(sigma_star_trans,1,z_LENGTH).*z-...
%         repmat(sigma_star_trans,1,z_LENGTH).^2/2).*exp(-z.^2/2).*...
%         normcdf(1./transpose(repmat(v,z_LENGTH,1)).*...
%         (z-repmat(sigma_star_trans./2,1,z_LENGTH)));
%         
% %     normcdf(1./v.*(z-repmat(sigma_star_trans./2,1,z_LENGTH)));
% expected_delta = 1/sqrt(2*pi)*sum(fz,2).*z_step;
% % For n=10 opt. FitGain
% n=10;
% c_n_10 = -n/2*transpose(log(1-2.*expected_delta./n))./p_eval;
% % For n=100 opt. FitGain
% n=100;
% c_n_100 = -n/2*transpose(log(1-2.*expected_delta./n))./p_eval;
% 

z_start = 0;
z_end = 20;
z_step = 10.^(-4);

v_LENGTH = length(exp_v_array);
c_n_10 = zeros(v_LENGTH,1);
c_n_100 = zeros(v_LENGTH,1);


for i = 1:1:v_LENGTH
    sigma_star = sigma_star_array(i);
    v = exp_v_array(i);
    z = (z_start:z_step:z_end)+sigma_star./2;
    p_eval = normcdf(-sigma_star./2./sqrt(1+v.^2));
    
    n=10;
    fz_10 = log(1-2.*((sigma_star*z-sigma_star/2).*exp(-z.^2/2).*normcdf(1/v*(z-sigma_star/2)))./n);
    E_delta_10 = 1/sqrt(2*pi)*sum(fz_10,2).*z_step;
    c_n_10(i) = -n/2*E_delta_10./p_eval;
    
    n=100;
    fz_100 = log(1-2.*((sigma_star*z-sigma_star/2).*exp(-z.^2/2).*normcdf(1/v*(z-sigma_star/2)))./n);
    E_delta_100 = 1/sqrt(2*pi)*sum(fz_100,2).*z_step;
    c_n_100(i) = -n/2*E_delta_100./p_eval;
    fprintf('%d,%d \n',c_n_10(i),c_n_100(i));
end

% for i = 1:1:v_LENGTH
%     sigma_star = sigma_star_array(i);
%     v = exp_v_array(i);
%     z = (z_start:z_step:z_end)+sigma_star./2;
%     p_eval = normcdf(-sigma_star./2./sqrt(1+v.^2));
%     fz = (sigma_star*z-sigma_star/2).*exp(-z.^2/2).*normcdf(1/v*(z-sigma_star/2));
%     E_delta = 1/sqrt(2*pi)*sum(fz,2).*z_step;
%     n=10;
%     c_n_10(i) = -n/2*log(1-2.*E_delta./n)./p_eval;
%     n=100;
%     c_n_100(i) = -n/2*log(1-2.*E_delta./n)./p_eval;
%     
% end

figure(1);
subplot(1,2,1);
scatterColour = 'r';

n=10;
d = sprintf('N=%d',n);
typeDot = 'x';
scatter(exp_v_array,c_n_10,typeDot,scatterColour,'DisplayName',d); hold on; 

n=100;
d = sprintf('N=%d',n);
typeDot = 'o';
scatter(exp_v_array,c_n_100,typeDot,scatterColour,'DisplayName',d); hold on;
set(gca, 'XScale', 'log');
