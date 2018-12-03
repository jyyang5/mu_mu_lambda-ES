% final version
% close all;
NUM_OF_RUNS = 50;

% oneLambda_RunALL(0,NUM_OF_RUNS);
FigureNum_opt1 = 10;    % Fig. numeber for first fitGain 
FigureNum_opt2 = 11;    % opt. sigmaStar and fitGain over v
lambda10 = oneLambda_RunALL_eta(1,FigureNum_opt1,FigureNum_opt2);
lambda20 = oneLambda_RunALL_eta(2,FigureNum_opt1,FigureNum_opt2);
lambda40 = oneLambda_RunALL_eta(3,FigureNum_opt1,FigureNum_opt2);

