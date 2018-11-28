% final version
% % close all;
% NUM_OF_RUNS = 50;
% 
% % oneLambda_RunALL(0,NUM_OF_RUNS);
% FigureNum_opt1 = 10;    % Fig. numeber for first fitGain 
% FigureNum_opt2 = 11;    % opt. sigmaStar and fitGain over v
% lambda10 = oneLambda_RunALL_eta(1,FigureNum_opt1,FigureNum_opt2);
% lambda20 = oneLambda_RunALL_eta(2,FigureNum_opt1,FigureNum_opt2);
% lambda40 = oneLambda_RunALL_eta(3,FigureNum_opt1,FigureNum_opt2);
% 



FigureNum = 1;

opt_fitGain_10 = 0.1703;
opt_fitGain_20 = 0.1844;
opt_fitGain_40 = 0.1929;
opt_fitGain_1_1 = 0.202;

opt_sigmastar_10 = 3.2001;
opt_sigmastar_20 = 6.0801;
opt_sigmastar_40 = 12.4201;
opt_sigmastar_1_1 = 1.224;

dotType = ':';
% Opt. FitGain
figure(FigureNum)
subplot(1,2,1);

legend('-DynamicLegend'); 
d1 = sprintf('N \\rightarrow \\infty');
figure(FigureNum);hold on;
subplot(1,2,1);
plot([0.1 10],[opt_fitGain_10 opt_fitGain_10],dotType,'Color','k','DisplayName',d1);hold on;
plot([0.1 10],[opt_fitGain_20 opt_fitGain_20],dotType,'Color','b','DisplayName',d1);hold on;
plot([0.1 10],[opt_fitGain_40 opt_fitGain_40],dotType,'Color','m','DisplayName',d1);hold on;
plot([0.1 10],[opt_fitGain_1_1 opt_fitGain_1_1],dotType,'Color','r','DisplayName',d1);hold on;

ylim([0 8]);   % set y = 0-10

legend('-DynamicLegend'); 
legend('show');
leg = legend();
leg.FontSize = 10;
title(leg,'Dimension of data');
set(gca,'FontSize',15);
box on;

legend('-DynamicLegend'); 
figure(FigureNum);hold on;
subplot(1,2,2);
d1 = sprintf('(3/3,10)-ES');
plot([0.1 10],[opt_sigmastar_10 opt_sigmastar_10],dotType,'Color','k','DisplayName',d1);hold on;
d1 = sprintf('(5/5,20)-ES');
plot([0.1 10],[opt_sigmastar_20 opt_sigmastar_20],dotType,'Color','b','DisplayName',d1);hold on;
d1 = sprintf('(10/10,40)-ES');
plot([0.1 10],[opt_sigmastar_40 opt_sigmastar_40],dotType,'Color','m','DisplayName',d1);hold on;
d1 = sprintf('(1+1)-ES');
plot([0.1 10],[opt_sigmastar_1_1 opt_sigmastar_1_1],dotType,'Color','r','DisplayName',d1);hold on;

ylim([0 15]);   % set y = 0-10
legend('-DynamicLegend'); 
legend('show');
leg = legend();
leg.FontSize = 10;
title(leg,'Strategy');
set(gca,'FontSize',15);
box on;
