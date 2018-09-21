% C_{mu/mu,lambda}
% for C_{3/3,10}=1.065389626877247
lambda = 10;
mu = 3;
step = 0.0000001;
x = -10:step:10;
int = exp(-x.^2).*(normcdf(x)).^(lambda-mu-1).*(1-normcdf(x)).^(mu-1);
% for i=-10:step:10
%     sum = sum + step*exp(-i^2)*(normcdf(i))^(lambda-mu-1)*(1-normcdf(i))^(mu-1);
% end
result = (lambda-mu)/(2*pi)*nchoosek(lambda,mu)*sum(int)*step;
disp(result);
