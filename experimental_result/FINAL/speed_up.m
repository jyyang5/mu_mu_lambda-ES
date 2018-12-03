

% Calculate speed-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% med_array = median(T_array,3);
% 
% median(T_array(),3)./repmat(T,1,4);
repmat(T,1,3)./median(T_array_8(:,2:4,:),3)

repmat(T,1,3)./median(T_array_16(:,2:4,:),3)

repmat(T,1,3)./median(T_array_32(:,2:4,:),3)


