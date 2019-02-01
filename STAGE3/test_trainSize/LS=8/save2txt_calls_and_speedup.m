% T_med_f6(:,[1,11],:)'
% T_med_f7(:,[1,11],:)'
% T_med_f8(:,[1,11],:)'
% T_med_f9(:,:)'


% fprintf('----------f6---------');
% temp_array6
% fprintf('----------f7---------');
% temp_array7
% fprintf('----------f8---------');
% temp_array8
% fprintf('----------f9---------');
% temp_array9


for j = 1:4
    fileID6 = fopen('f6_data.txt','a+');
    fileID7 = fopen('f7_data.txt','a+');
    fileID8 = fopen('f8_data.txt','a+');
    
    % Dim esperimented
    n_array = [2,4,8,16];
    % Index arrays want to stored for each fun
    p_f6_array = [1, 6, 9, 10, 11];
    p_f7_array = [1, 11];
    p_f8_array = [2,10];

    
    n = n_array(j);
    if n==2
        load('all_dim=2.mat');
    elseif n==4
        load('all_dim=4.mat')
    elseif n==8
        load('all_dim=8.mat')
    elseif n==16
        load('all_dim=16.mat')
    end
    % Stores T 
    temp_array6 = zeros(length(p_f6_array),9);
    temp_array7 = zeros(length(p_f7_array),9);
    temp_array8 = zeros(length(p_f8_array),9);


    temp_array6(:,[1,2,4,6,8]) = T_med_f6(:,p_f6_array,:)';
    temp_array6(:,3:2:9) = repmat(temp_array6(:,1),1,4)./temp_array6(:,[2,4,6,8]);

    temp_array7(:,[1,2,4,6,8]) = T_med_f7(:,p_f7_array,:)';
    temp_array7(:,3:2:9) = repmat(temp_array7(:,1),1,4)./temp_array7(:,[2,4,6,8]);

    temp_array8(:,[1,2,4,6,8]) = T_med_f8(:,p_f8_array,:)';
    temp_array8(:,3:2:9) = repmat(temp_array8(:,1),1,4)./temp_array8(:,[2,4,6,8]);


    % Write file for f6[sphere function]
    for p = 1:length(p_f6_array)
        para_index = p_f6_array(p);
        for i = 1:5
            if i==1
                fprintf(fileID6, sprintf('&    %d &    %.2f &    %d &    ',n,f6_range(para_index),temp_array6(p,i)));
            elseif i==5
                fprintf(fileID6, '%d (%.2f) \\\\ \n    ',temp_array6(p,2*i-2),temp_array6(p,2*i-1));
            else
                fprintf(fileID6, '%d (%.2f) & %    ',temp_array6(p,2*i-2),temp_array6(p,2*i-1));
            end
        end
        fprintf(fileID6, '\n');
    end
    % Write file for f7[quartic function]
    for p = 1:length(p_f7_array)
        para_index = p_f7_array(p);
        for i = 1:5
            if i==1
                fprintf(fileID7, sprintf('&    %d &    %.2f &    %d &    ',n,f7_range(para_index),temp_array7(p,i)));
            elseif i==5
                fprintf(fileID7, '%d (%.2f) \\\\ \n    ',temp_array7(p,2*i-2),temp_array7(p,2*i-1));
            else
                fprintf(fileID7, '%d (%.2f) & %    ',temp_array7(p,2*i-2),temp_array7(p,2*i-1));
            end
        end
        fprintf(fileID7, '\n');
    end
    % Write file for f8[elliposids function]
    for p = 1:length(p_f8_array)
        para_index = p_f8_array(p);
        for i = 1:5
            if i==1
                fprintf(fileID8, sprintf('&    %d &    %.2f &    %d &    ',n,f8_range(para_index),temp_array8(p,i)));
            elseif i==5
                fprintf(fileID8, '%d (%.2f) \\\\ \n    ',temp_array8(p,2*i-2),temp_array8(p,2*i-1));
            else
                fprintf(fileID8, '%d (%.2f) & %    ',temp_array8(p,2*i-2),temp_array8(p,2*i-1));
            end
        end
        fprintf(fileID8, '\n');
    end
    
end


fclose(fileID6);
fclose(fileID7);
fclose(fileID8);
% fclose(fileID9);






% 
% for i = 1:5
%     if i==1
%         temp_array6(1,1) = T_med_f6(1,1,:)';
%         temp_array6(2,1) = T_med_f6(1,11,:)';
%         temp_array7(1,1) = T_med_f7(1,1,:)';
%         temp_array7(2,1) = T_med_f7(1,11,:)';
%         temp_array8(1,1) = T_med_f8(1,1,:)';
%         temp_array8(2,1) = T_med_f8(1,11,:)';
%     else
%         temp_array6(1,2*i) = T_med_f6(i,1,:)';
%         temp_array6(1,2*i+1) = T_med_f6(i,1,:)'/;
%         temp_array6(2,2*i) = T_med_f6(i,11,:)';
%         temp_array7(1,2*i+1) = T_med_f7(i,1,:)';
%         temp_array7(2,2*i) = T_med_f7(i,11,:)';
%         temp_array8(1,2*i+1) = T_med_f8(i,1,:)';
%         temp_array8(2,i) = T_med_f8(i,11,:)';
%     end
%     
%     
% end