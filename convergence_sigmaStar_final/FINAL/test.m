% final version
close all;
NUM_OF_RUNS = 10;
f = @(x) (x'*x);
FIG_NUM = 2;
subplotNum_step = 2;
% s_start = 0.0001;
% increment =1;
% s_end = 10+s_start;
% v_array = s_start:increment:s_end;                                          % For experimental result
opt_plot_colour = 'r';
c_mu_lambda = 0;
v_expedted_curve_array = exp(-2.302585092994046: 0.0461:2.302585092994046+0.01);
v_array = v_expedted_curve_array(1:10:101);

fun_precise_optFitGain_over_v_ONE(f,NUM_OF_RUNS,FIG_NUM,subplotNum_step,v_array,10,opt_plot_colour,'x',c_mu_lambda,v_expedted_curve_array);
fun_precise_optFitGain_over_v_ONE(f,NUM_OF_RUNS,FIG_NUM,subplotNum_step,v_array,100,opt_plot_colour,'o',c_mu_lambda,v_expedted_curve_array);