n_init = lambda;
nb = max(1, floor(lambda/10));
offspring = zeros(n,mu);			% mu offspring chosen
%candidate_unevaluated = zeros(n,lambda);
fep_offspring = zeros(mu,1);
fep_temp = zeros(lambda,1);
xTrain = zeros(10,10000);
fTrain = zeros(1,10000);
T = 0;						% number of objective function calls

% evaluate use GP
for i=1:1:lambda
		z(:,i) = randn(n,1);
        y(:,i) = centroid + sigma*z(:,i);
        fyep(i) = gp(xTrain, fTrain, y(:,i), theta);
end
% for simple calculation 
fy = fyep;
% sort fyep (smaller first)
[index, initial_sorted_order] = sort(fy);
y = y(:,initial_sorted_order);
% update training set
offspring = y(:,1:mu);
xTrain(:,T:T+mu-1) = offspring;
for i = 1:1:mu
	fTrain(T) = f(offspring(:,i));
	T = T + 1;
end
% set for unevaluated candidate solutions
candidate_unevaluated =  y(:,mu+1:lambda);

for i=1:1:(lambda-n_init)/nb
	
	% Reevaluate offsprings
	for j=1:1:mu
		fep_offspring(j) = gp(xTrain(:,T-40:T-1), fTrain(T-40:T-1), offspring(:,j), theta);
	end
	[index, sorted_fep] = sort(fep_offspring);
	% order of mu candidate solutions unchanged
	if (sorted_fep(1:mu) == 1:1:mu)
		break;
	% ranking change (i.e. GP not stable)
	else
		% add mu bset unevaluated points to GP training set
		for j=1:1:lambda-mu*i
			fep_temp(j) = gp(xTrain(:,T-40:T-1), fTrain(T-40:T-1), candidate_unevaluated(:,j), theta);
		end
		fep_temp = fep_temp(1:lambda-mu*i);					% change size of the fep array 
		[index, sorted_fep] = sort(fep_temp);
		candidate_unevaluated = candidate_unevaluated(:,sorted_fep);
		% add next nb unevaluated points to GP training 
		offspring = candidate_unevaluated(:,1:nb);
		xTrain(:,T+nb-1) = offspring;
		for j=1:1:nb
			fTrain(T) = f(xTrain(:,T));
			T = T + 1;
		end
		
		[dim, len_unevaluated] = size(candidate_unevaluated); 
		candidate_unevaluated = candidate_unevaluated(:,nb+1:len_unevaluated); 
	end
	if (i>2)
		n_init = min(n_init,lambda-nb);
	elseif (i<2)
		n_init = max(nb,n_init-nb);
end 
