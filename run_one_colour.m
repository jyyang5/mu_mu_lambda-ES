% run fun_fitness_sigmaStar to plot convergence rate 
% at noise-to-signal ratio 0,.25,1,4 respectively
% also for data with dim = 10, 100

%fun_fitness_sigmaStar(ita,n,scatterColour,typeDot,LINE_OR_NOT)
%Input
%   mu                    # parents AND # of parents replaced by offspring
%   lambda                generated offspring size
%   ita:                  ita = sigma_ep_star /sigma_star
%   n:                    dim of data
%   scatterColour:        r  Red, g  Green, b   Blue, c   Cyan, m  Magenta, y  Yellow, k   Black,w  White
%   typeDot:              type of dot (+  Plus sign, o  Circle, *  Asterisk, x    Cross)
%   LINE_OR_NOT:          if just scatter = 0
%Return
%   scatter or scatter+curve

%
mu = 3;
lambda = 10;

figure(4);
% green
ita = 2;
%fun_fitness_sigmaStar(4,100000,'g','.',1);        
fun_fitness_sigmaStar(mu,lambda,ita,10,'g','x',0);        
fun_fitness_sigmaStar(mu,lambda,ita,100,'g','o',0);

disp('green done');


% red
ita = 1;
fun_fitness_sigmaStar(mu,lambda,ita,10,'r','x',0);        
fun_fitness_sigmaStar(mu,lambda,ita,100,'r','o',0);

disp('red done');
% blue
ita = 0.25;
fun_fitness_sigmaStar(mu,lambda,ita,10,'b','x',0);        
fun_fitness_sigmaStar(mu,lambda,ita,100,'b','o',0);

disp('blue done');


% black
ita = 0;
fun_fitness_sigmaStar(mu,lambda,ita,10,'k','x',0);        
fun_fitness_sigmaStar(mu,lambda,ita,100,'k','o',0);

disp('black done');