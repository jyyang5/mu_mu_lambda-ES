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
    fileID9 = fopen('f9_data.txt','a+');

    n_array = [2,4,8,16];
    p_array = [1, 11];

    
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
    temp_array6 = zeros(2,9);
    temp_array7 = zeros(2,9);
    temp_array8 = zeros(2,9);
    temp_array9 = zeros(1,9);


    temp_array6(:,[1,2,4,6,8]) = T_med_f6(:,p_array,:)';
    temp_array6(:,3:2:9) = repmat(temp_array6(:,1),1,4)./temp_array6(:,[2,4,6,8]);

    temp_array7(:,[1,2,4,6,8]) = T_med_f7(:,p_array,:)';
    temp_array7(:,3:2:9) = repmat(temp_array7(:,1),1,4)./temp_array7(:,[2,4,6,8]);

    temp_array8(:,[1,2,4,6,8]) = T_med_f8(:,p_array,:)';
    temp_array8(:,3:2:9) = repmat(temp_array8(:,1),1,4)./temp_array8(:,[2,4,6,8]);

    temp_array9(:,[1,2,4,6,8]) = T_med_f9(:,:)';
    temp_array9(:,3:2:9) = repmat(temp_array9(:,1),1,4)./temp_array9(:,[2,4,6,8]);

    for p = 1:length(p_array)
        para_index = p_array(p);
        for i = 1:5
            if i==1

                fprintf(fileID6, sprintf('&    %d &    %.2f &    %d &    ',n,f6_range(para_index),temp_array6(p,i)));
                fprintf(fileID7, sprintf('&    %d &    %.2f &    %d &    ',n,f7_range(para_index),temp_array7(p,i)));
                fprintf(fileID8, sprintf('&    %d &    %.2f &    %d &    ',n,f8_range(para_index),temp_array8(p,i)));
                if p==1
                    fprintf(fileID9, sprintf('&    %d &     &    %d &    ',n,temp_array9(p,i)));
                end
            elseif i==5
                fprintf(fileID6, '%d (%.2f) \\\\ \n    ',temp_array6(p,2*i-2),temp_array6(p,2*i-1));
                fprintf(fileID7, '%d (%.2f) \\\\ \n    ',temp_array7(p,2*i-2),temp_array7(p,2*i-1));
                fprintf(fileID8, '%d (%.2f) \\\\ \n    ',temp_array8(p,2*i-2),temp_array8(p,2*i-1));
                if p==1
                    fprintf(fileID9, '%d (%.2f) \\\\ \n    ',temp_array9(p,2*i-2),temp_array9(p,2*i-1));
                end
            else
                fprintf(fileID6, '%d (%.2f) & %    ',temp_array6(p,2*i-2),temp_array6(p,2*i-1));
                fprintf(fileID7, '%d (%.2f) & %    ',temp_array7(p,2*i-2),temp_array7(p,2*i-1));
                fprintf(fileID8, '%d (%.2f) & %    ',temp_array8(p,2*i-2),temp_array8(p,2*i-1));
                if p==1
                    fprintf(fileID9, '%d (%.2f) & %    ',temp_array9(p,2*i-2),temp_array9(p,2*i-1));
                end
            end
        end
        fprintf(fileID6, '\n');
        fprintf(fileID7, '\n');
        fprintf(fileID8, '\n');
        if p==1
            fprintf(fileID9, '\n');
        end
    end
end

fclose(fileID6);
fclose(fileID7);
fclose(fileID8);
fclose(fileID9);






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