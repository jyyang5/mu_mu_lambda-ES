% Plot solid line in Expected fitness gain over sigma* and sigma_ep^* 
% solid line in convergence rate graph

a = openfig('test.fig');

mu = 3;
s_start = 0.2;
increment =0.2;
s_end = 6.2;
sigma_star = s_start:increment:s_end;

% noise-to-signal ratio = 4
ita = 4;
final = sigma_star*(1.065389626877247)./sqrt(1+ita^2)-sigma_star.^2./(2*mu);
% noise-to-signal ratio = 1
ita = 1;
final1 = sigma_star*(1.065389626877247)./sqrt(1+ita^2)-sigma_star.^2./(2*mu);
% noise-to-signal ratio = 1/4
ita = 1/4;
final2 = sigma_star*(1.065389626877247)./sqrt(1+ita^2)-sigma_star.^2./(2*mu);
% noise-to-signal ratio = 0
ita = 0;
final3 = sigma_star*(1.065389626877247)./sqrt(1+ita^2)-sigma_star.^2./(2*mu);

figure(2);
hold on;
plot(sigma_star,final,'g');     % ita = 4
plot(sigma_star,final1,'r');    % ita = 1
plot(sigma_star,final2,'b');    % ita = 1/4
plot(sigma_star,final3,'k');    % ita = 0
hold off;




