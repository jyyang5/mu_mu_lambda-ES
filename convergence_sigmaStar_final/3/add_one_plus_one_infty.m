% First load data .mat
legend('-DynamicLegend');
figure(2);
subplot(1,2,1);hold on; 
d = sprintf('N \\rightarrow \\infty');
scatterColour = 'r';
plot(v_array,opt_eta_array,scatterColour,'DisplayName',d);
legend('-DynamicLegend'); 
legend('show');
leg = legend();
leg.FontSize = 10;


legend('-DynamicLegend');
figure(2);
subplot(1,2,2);hold on; 
d = sprintf('N \\rightarrow \\infty');
scatterColour = 'r';
plot(v_array,opt_sigma_star_array,scatterColour,'DisplayName',d);
legend('-DynamicLegend'); 
legend('show');
leg = legend();
leg.FontSize = 10;
