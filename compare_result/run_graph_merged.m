% test file
% using 5 different obejctive functions plot 5 seperate graphs
% NOTE: DO NOT CHANGE THE NUMBER WHEN CALL fun_graph_merged
%       OR CANNOT SAVE ITERTAION DATA

f1 = @(x) (x'*x)^(1/2);
f2 = @(x) (x'*x);
f3 = @(x) (x'*x)^(3/2);


NUM_OF_RUNS = 400;


disp("===========================================================");
disp("linear sphere");
disp("---------------");
%graph_fun(f1,1);
T_f1=fun_graph_merged(f1,6,NUM_OF_RUNS);

disp("===========================================================");
disp("quadratic sphere");
disp("---------------");
T_f2=fun_graph_merged(f2,7,NUM_OF_RUNS);

disp("===========================================================");
disp("cubic sphere");
disp("---------------");
T_f3=fun_graph_merged(f3,8,NUM_OF_RUNS);

disp("===========================================================");
disp("schwefel's function");
disp("---------------");
T_f4=fun_graph_merged(@f4,9,NUM_OF_RUNS);

disp("===========================================================");
disp("quartic function");
disp("---------------");
T_f5=fun_graph_merged(@f5,10,NUM_OF_RUNS);
disp("===========================================================");



function [val] = f4(x)

val = 0;
for i = 1:length(x)
    temp = 0;
    for j = 1:i
        temp = temp + x(j);
    end
    val = val + temp^2;    
        
end


end

function [val] = f5(x)

val = 0;
beta = 1;
for i = 1:length(x)-1
    val = val + beta*(x(i+1)-x(i)^2)^2+(1-x(i))^2;
end
end