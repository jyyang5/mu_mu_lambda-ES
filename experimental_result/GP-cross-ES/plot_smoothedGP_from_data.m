


TRAINING_SIZE = 40;

  
% GP smooth 
window_length = 40;
kernel = exp(-(-3*window_length:3*window_length).^2/window_length^2/2);
kernel = kernel/sum(kernel);        % Normalized 

sigma = 40;
kernel = exp(-[-3*sigma:3*sigma].^2/sigma^2/2);
kernel = kernel/sum(kernel);

for fname = 1:1:5
    for i = 1:1:length(lambda_array)
        % Find the index of median T 
        T_array_temp = T_array(fname,i,:);
        sorted_T = sort(T_array_temp);
        med_index = find(T_array_temp == sorted_T(ceil(length(sorted_T)/2)));
        med_index = med_index(1);   % randomly pick one
        T_med = T_array(fname,i,med_index);
        GP_error_matrix_med = transpose(squeeze(GP_error_matrix(fname,i,med_index,:)));
        t_med = squeeze(t_array(fname,i,med_index));
        f_x_med = transpose(squeeze((f_x_matrix(fname,i,med_index,:))));
        
        lambda =  lambda_array(i);
        figure(300);
        legend('-DynamicLegend'); 
        hold on;
        mu = ceil(lambda/4);
        if lambda==0
            d = sprintf('(1+1)-ES');
        else
            d = sprintf('(%d/%d,%d)-ES',mu,mu,lambda);
        end
        subplot(1,subplot_COL,fname);
        
        % GP smooth 
        window_length = 40;
        kernel = exp(-(-3*window_length:3*window_length).^2/window_length^2/2);
        kernel = kernel/sum(kernel);        % Normalized 
        
        if lambda==0
            t_start = TRAINING_SIZE+3;
            plot(t_start:1:T_med,GP_error_matrix_med(t_start:T_med),'DisplayName',d);hold on;
%             smoothed_GP = smoothdata(GP_error_matrix_med(t_start:T_med),'gaussian',40);
                
            smoothed_GP = exp(conv(log(GP_error_matrix_med(t_start:t_med)), kernel, 'same'));
            d1 = sprintf('(1+1)-ES[smoothed]');
            plot(t_start:1:t_med,smoothed_GP,'DisplayName',d1,'LineWidth',2);hold on;
%             set(gca, 'YScale', 'log');
        else 
%             GP_error_range1 = GP_error_matrix_med(1:t_start);
            t_start = ceil(TRAINING_SIZE/lambda);
            fx_range2 = f_x_med(t_start+1:t_med);
            t_range2 = lambda*t_start+1:lambda*t_start+length(fx_range2);
            GP_error_range2 = GP_error_matrix_med(t_start+1:t_med);
            plot(t_range2,GP_error_range2,'DisplayName',d);hold on;
%             plot([t_range1 t_range2], [GP_error_range1 GP_error_range2],'DisplayName',d);hold on;
            % GP smoothed
%             smoothed_GP_range1 = smoothdata(GP_error_range1,'gaussian',40);
%             smoothed_GP_range2 = smoothdata(GP_error_range2,'gaussian',40);
            smoothed_GP_range2 = exp(conv(log(GP_error_range2), kernel, 'same'));

            d1 = sprintf('(%d/%d,%d)-ES[smoothed]',mu,mu,lambda);
            plot(t_range2, smoothed_GP_range2,'DisplayName',d1,'LineWidth',2);hold on;
            set(gca, 'YScale', 'log');

%             plot([t_range1 t_range2],[smoothed_GP_range1 smoothed_GP_range2],'DisplayName',d1,'LineWidth',2);hold on;
        end
%         ylim([10^(-2) 10^4]);

    end
end
