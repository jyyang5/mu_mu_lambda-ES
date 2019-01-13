PROB_RATE_array = 0.1:0.1:0.5;
legendCell = {};
for i = 1:1:length(PROB_RATE_array)
    legendCell{i} = sprintf('C3=%.2f',PROB_RATE_array(i));
%     histogram(randn(10,1));hold on;
end

legend(legendCell);

% legendCell = {};
% for i = 1:1:length(PROB_RATE_array)
%     legendCell(i) = sprintf('C3=%.2f',PROB_RATE_array(i));
% end
% legend(legendCell);