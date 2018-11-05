% test file
% using 5 different obejctive functions plot 5 seperate graphs
% NOTE: DO NOT CHANGE THE NUMBER WHEN CALL fun_graph_merged
%       OR CANNOT SAVE ITERTAION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = @(x) (x'*x)^(1/2);
f2 = @(x) (x'*x);
f3 = @(x) (x'*x)^(3/2);

close all;
% NUM_OF_RUNS = 1000;
% % NUM_OF_RUNS = 2;
% TRAINING_SIZE = 40;
% LENGTH_SCALE = 8;

subplot_ROW = 5;
subplot_COL = 5; 
% 
% t_array = zeros(5,4,NUM_OF_RUNS,1);                                           % # of iterations for the stop creteria
% sigma_matrix = zeros(5,4,NUM_OF_RUNS,10000);                                  % store all sigma
% T_array = zeros(5,4,NUM_OF_RUNS,1);                                           % # of objective function evaluations for the stop creteria
% f_x_matrix = zeros(5,4,NUM_OF_RUNS,10000);                                    % store all fx
% convergence_rate_array = zeros(5,4,NUM_OF_RUNS,1);                            % convergence rate 
% GP_error_matrix = zeros(5,4,NUM_OF_RUNS,10000);                               % store similar to noise-to-signal ratio
% sigma_star_matrix = zeros(5,4,NUM_OF_RUNS,10000);                             % normalized step size 
% success_rate_array = zeros(5,4,NUM_OF_RUNS,1);                                % success rate 
% delta_matrix = zeros(5,4,NUM_OF_RUNS,10000);                                  % each [i,j] stores a delta array 

% t_med = median(t_array);
% t_max = max(t_array);
% t_min = min(t_array);
% sigma_matrix_med = squeeze(median(sigma_matrix));
% T_med = median(T_array);
% f_x_med = squeeze(median(f_x_matrix));
% convergence_med = median(convergence_rate_array);
% GP_error_matrix_med = squeeze(median(GP_error_matrix));
% sigma_star_matrix_med = squeeze(median(sigma_star_matrix));
% success_med = median(success_rate_array);
% delta_matrix_med = squeeze(median(delta_matrix));
lambda_array = [0,10,20,40];

close all;
for fname = 1:1:5
    for i = 1:1:length(lambda_array)
        
        t_med = median(t_array(fname,i,:));
        t_max = max(t_array(fname,i,:));
        t_min = min(t_array(fname,i,:));
        sigma_matrix_med = median(squeeze(sigma_matrix(fname,i,:,:)));
        T_med = median(T_array(fname,i,:));
        f_x_med = median(squeeze((f_x_matrix(fname,i,:,:))));
        convergence_med = median(squeeze((convergence_rate_array(fname,i,:))));
        GP_error_matrix_med = median(squeeze(GP_error_matrix(fname,i,:,:)));
        sigma_star_matrix_med = median(squeeze(sigma_star_matrix(fname,i,:,:)));
        success_med = median(success_rate_array(fname,i,:));
        delta_matrix_med = median(squeeze(delta_matrix(fname,i,:,:)));
        
        figure(200);
        legend('-DynamicLegend'); 
        hold on;
        mu = ceil(lambda/4);
        if lambda==0
            d = sprintf('(1+1)-ES');
        else
            d = sprintf('(%d/%d,%d)-ES',mu,mu,lambda);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1.objective function [row 1]
        subplot(subplot_ROW,subplot_COL,fname);
        h = histogram(T_array,'Normalization','probability','DisplayName',d);hold on;
%         h.BinWidth = 20;          % Bin width
        if(lambda==40)
            legend({'(1+1)-ES','(3/3,10)-ES','(5/5,20)-ES','(10/10,40)-ES'})
            if(fname == 1)
                d3 =sprintf('linear sphere');
                title(d3,'fontsize',15);
            elseif(fname == 2)
                d3 =sprintf('quadratic sphere');
                title(d3,'fontsize',15);
            elseif(fname == 3)
                d3 =sprintf('cubic sphere');
                title(d3,'fontsize',15);
            elseif(fname == 4)
                d3 =sprintf('Schwefel;s function');
                title(d3,'fontsize',15);
            elseif(fname == 5)
                d3 =sprintf('quartic function');
                title(d3,'fontsize',15); 
            end
        end
        xlabel('objective function calls','FontSize',15); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2.objective function [row 2]
        subplot(subplot_ROW,subplot_COL,fname+5);
        if lambda==0
            t_start = TRAINING_SIZE+1+2;
            plot(1:1:t_med,f_x_med(1:T_med),'DisplayName',d);hold on;
        else 
            t_start = ceil(TRAINING_SIZE/lambda);
            fx_range1 = f_x_med(1:t_start);
            fx_range2 = f_x_med(t_start+1:t_med);
            t_range1 = 1:lambda:lambda*t_start;
            t_range2 = lambda*t_start+1:lambda*t_start+length(fx_range2);
            plot([t_range1 t_range2], [fx_range1 fx_range2],'DisplayName',d);hold on;% mml with GP
        end
        if(fname==1)
            ylabel('objective function value','FontSize',15);%
        end
        xlabel('objective function calls','FontSize',15); 
        set(gca, 'YScale', 'log');
        legend('-DynamicLegend'); 
        legend('show');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3.GP error [row 3]
        subplot(subplot_ROW,subplot_COL,fname+10);
        if lambda==0
            plot(t_start:1:T_med,GP_error_matrix_med(t_start:T_med),'DisplayName',d);hold on;
            d1 = sprintf('(%d/%d,%d)-ES[smoothed]',mu,mu,lambda);
            plot(t_start:1:T_med,smoothdata(GP_error_matrix_med(t_start:T_med),'gaussian',40),'DisplayName',d1,'LineWidth',2);hold on;
        else 
            GP_error_range1 = GP_error_matrix_med(1:t_start);
            GP_error_range2 = GP_error_matrix_med(t_start+1:t_med);
            plot([t_range1 t_range2], [GP_error_range1 GP_error_range2],'DisplayName',d);hold on;
            % GP smoothed
            smoothed_GP_range1 = smoothdata(GP_error_range1,'gaussian',40);
            smoothed_GP_range2 = smoothdata(GP_error_range2,'gaussian',40);
            d1 = sprintf('(%d/%d,%d)-ES[smoothed]',mu,mu,lambda);
            plot([t_range1 t_range2],[smoothed_GP_range1 smoothed_GP_range2],'DisplayName',d1,'LineWidth',2);hold on;
        end
        if(fname==1)
            ylabel('relative model error','FontSize',15);%
        end
        xlabel('objective function evaluations','FontSize',15); 
        set(gca, 'YScale', 'log');

        legend('-DynamicLegend'); 
        legend('show');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4.normalized step size [row 4](only makes sense for sphere functions)
        if(fname<4)
            subplot(subplot_ROW,subplot_COL,fname+15);
            if lambda==0
                plot(1:1:T_med,sigma_star_matrix_med(1:T_med),'DisplayName',d);hold on;
            else 
                sigma_star_range1 = sigma_star_matrix_med(1:t_start);
                sigma_star_range2 = sigma_star_matrix_med(t_start+1:t_med);
                plot([t_range1 t_range2], [sigma_star_range1 sigma_star_range2],'DisplayName',d);hold on;
            end
            if(fname==1)
                ylabel('normalized step size \sigma*','FontSize',15);%
            end
            xlabel('objective function evaluations','FontSize',15); 
            set(gca, 'YScale', 'log');
        %     title('normalized step size \sigma*','FontSize',20);

            legend('-DynamicLegend'); 
            legend('show');
        end
        % 5. histogram success rate (good step size)
        figure(300);
        % No bar colour, bold binwidth
        if(T_med == 5000) % early stopping
            h2 = histogram(nonzeros(delta_matrix_med),'Normalization','probability','DisplayName',d);hold on;
%             h2 = histogram(nonzeros(delta_matrix_med),'Normalization','probability','DisplayName',d,'FaceColor','none','LineWidth',2);hold on;

        else
            h2 = histogram(nonzeros(delta_matrix_med),'Normalization','probability','DisplayName',d);hold on;
%             h2 = histogram(nonzeros(delta_matrix_med),'Normalization','probability','DisplayName',d,'FaceColor','none','LineWidth',2);hold on;
        end
        if(lambda==0)
            h2.EdgeColor= [0  0.4470 0.7410];
%             h2.BinEdges=h2.BinEdges-0.015;
        elseif(lambda == 10)
            h2.EdgeColor= [0.8500  0.3250  0.0980];
%             h2.BinEdges=h2.BinEdges-0.005;
        elseif(lambda==20)
            h2.EdgeColor= [0.9290  0.6940  0.1250];
%             h2.BinEdges=h2.BinEdges+0.005;
        elseif(lambda==40)
            h2.EdgeColor= [0.4940  0.1840  0.5560];
%             h2.BinEdges=h2.BinEdges+0.015;
        end
%         h2.LineWidth=1.2;
%         h2.BinWidth = 0.2;

        figure(200);
        subplot(subplot_ROW,subplot_COL,fname+20);
        value = h2.Values;		% height of the bar
        width = h2.BinWidth;				% width of the bar
        range = h2.BinLimits;		% [startX endX]
        plot(range(1)+width/2:width:range(2)-width/2,value);hold on;
        if(fname==1)
            ylabel('Probability','FontSize',15);%
        end
        xlabel('Normalized fitness gain','FontSize',15); 
        % set(gca, 'YScale', 'log');
        % title('step size \sigma','FontSize',20);
        legend('-DynamicLegend'); 
        legend('show');

        saveas(gcf,'merged_plot_v1.fig');  
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Success rate good step size
        figure(201);
        % 5. histogram success rate (good step size)
        subplot(1,subplot_COL,fname);
        % No bar colour, bold binwidth
        if(T_med == 5000) % early stopping
            h2 = histogram(success_rate_array,'Normalization','probability','DisplayName',d);hold on;
        else
            h2 = histogram(success_rate_array,'Normalization','probability','DisplayName',d);hold on;
        end
        % if(lambda==0)
        %     h2.BinEdges=h.BinEdges-0.0125;
        % elseif(lambda == 10)
        %     h2.BinEdges=h.BinEdges;
        % elseif(lambda==20)
        %     h2.BinEdges=h.BinEdges+0.0125;
        % elseif(lambda==40)
        %     h2.BinEdges=h.BinEdges+0.025;
        % end
        % h2.LineWidth=1.5;
        h2.BinWidth = 0.05;
        if(fname==1)
            ylabel('Probability','FontSize',15);%
        end
        xlabel('Prob good step size','FontSize',15); 
        % set(gca, 'YScale', 'log');
        % title('step size \sigma','FontSize',20);
        legend('-DynamicLegend'); 
        legend('show');
        saveas(gcf,'prob_good_stepSize_v1.fig');
        

        
    end   
end

% save('data.mat','NUM_OF_RUNS','TRAINING_SIZE','LENGTH_SCALE','t_array','sigma_matrix',...
%     'T_array','f_x_matrix','convergence_rate_array','GP_error_matrix','sigma_star_matrix',...
%     'success_rate_array','delta_matrix');




% 
% %Get data
% a = cell2mat(temp1);
% b = cell2mat(temp2);
% c = cell2mat(temp3);
% % d = cell2mat(temp4);
% % e = cell2mat(temp5);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot summary
% figure(100)
% 
% % plot objective function evaluation
% subplot(3,1,3)
% hold on;
% plot(lambda_array,a(:,1));hold on;
% plot(lambda_array,b(:,1));hold on;
% plot(lambda_array,c(:,1));hold on;
% % plot(lambda_array,d(1));hold on;
% % plot(lambda_array,e(1));hold on;
% xlabel('Population size \lambda','fontsize',15);
% ylabel('Objective function evaluation','fontsize',15);
% legend('Linear sphere','Quadratic sphere','Cubic sphere');
% % legend('Linear sphere','Quadratic sphere','Cubic sphere','Schwefel?s function','Quartic function')
% d =sprintf('Objective function evaluation');
% title(d,'fontsize',20);
% 
% % plot convergence rate
% subplot(3,1,1)
% hold on;
% plot(lambda_array,a(:,2));hold on;
% plot(lambda_array,b(:,2));hold on;
% plot(lambda_array,c(:,2));hold on;
% % plot(lambda_array,d(2));hold on;
% % plot(lambda_array,e(2));hold on;
% xlabel('Population size \lambda','fontsize',15);
% ylabel('Convergence rate','fontsize',15);
% legend('Linear sphere','Quadratic sphere','Cubic sphere');%,'Schwefel?s function','Quartic function')
% d =sprintf('Convergence rate');
% title(d,'fontsize',20);
% 
% % plot success rate
% subplot(3,1,2)
% hold on;
% plot(lambda_array,a(:,3));hold on;
% plot(lambda_array,b(:,3));hold on;
% plot(lambda_array,c(:,3));hold on;
% % plot(lambda_array,d(3));hold on;
% % plot(lambda_array,e(3));hold on;
% xlabel('Population size \lambda','fontsize',15);
% ylabel('Objective function evaluation','fontsize',15);
% legend('Linear sphere','Quadratic sphere','Cubic sphere');
% % legend('Linear sphere','Quadratic sphere','Cubic sphere','Schwefel?s function','Quartic function')
% d =sprintf('Objective function evaluation');
% title(d,'fontsize',20);
% 
% % Saves fig
% saves(gcf,'summary.fig');
% % 
% % % linear sphere
% % disp("===========================================================");
% % disp("linear sphere");
% % disp("---------------");
% % name = 6;
% % result_f1=fun_graph_funCall_NEW(f1,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);
% % 
% % % quadratic sphere
% % disp("===========================================================");
% % disp("quadratic sphere");
% % disp("---------------");
% % name = 7;
% % result_f2 = fun_graph_funCall_NEW(f2,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);
% % 
% % % cubic sphere
% % disp("===========================================================");
% % disp("cubic sphere");
% % disp("---------------");
% % name = 8;
% % result_f3=fun_graph_funCall_NEW(f3,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);
% 
% % disp("===========================================================");
% % disp("schwefel's function");
% % disp("---------------");
% % name = 9;
% % result_f4=fun_graph_funCall_NEW(@f4,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);
% % 
% % disp("===========================================================");
% % disp("quartic function");
% % disp("---------------");
% % name = 10;
% % result_f5=fun_graph_funCall_NEW(@f5,name,NUM_OF_RUNS,mu,lambda,TRAINING_SIZE);
% % disp("===========================================================");
% 
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Display statics
% % plotValue(result_f1,6);
% % plotValue(result_f2,7);
% % plotValue(result_f3,8);
% % plotValue(result_f4,9);
% % plotValue(result_f5,10);
% % % save('final_result.mat','result_f1','result_f2','result_f3','NUM_OF_RUNS','TRAINING_SIZE','mu','lambda');
% % save('final_result.mat','result_f1','result_f2','result_f3','result_f4','result_f5','NUM_OF_RUNS','TRAINING_SIZE','mu','lambda');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [val] = f4(x)
%     val = 0;
%     for i = 1:length(x)
%         temp = 0;
%         for j = 1:i
%             temp = temp + x(j);
%         end
%         val = val + temp^2;          
%     end
% end
% 
% function [val] = f5(x)
%     val = 0;
%     beta = 1;
%     for i = 1:length(x)-1
%     val = val + beta*(x(i+1)-x(i).^2)^2+(1-x(i))^2;
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function plotValue(temp,name)
%     if name == 6
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%         disp('Linear sphere');
%     elseif name == 7
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%         disp('Quadratic sphere');    
%     elseif name == 8
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%         disp('Cubic sphere'); 
%     elseif name == 9 
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%         disp('Schwefel function');    
%     elseif name == 10
%         disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%         disp('Quartic function');
%     end
%         
%         
%     T = cell2mat(temp(1));
%     T1 = cell2mat(temp(2));
%     convergence = cell2mat(temp(3));
%     convergence1 = cell2mat(temp(4));
%     success = cell2mat(temp(5));
%     success1 = cell2mat(temp(6));
%     disp('# objective function evaluations (with model assistance)');
%     d = sprintf('mml-ES: %d',T);
%     disp(d);
%     d = sprintf('(1+1)-ES: %d',T1);
%     disp(d);
%     disp('Convergence rate (with model assistance)');
%     d = sprintf('mml-ES: %.3f',convergence);
%     disp(d);
%     d = sprintf('(1+1)-ES: %.3f',convergence1);
%     disp(d);
%     disp('Success rate (with model assistance)');
%     d = sprintf('mml-ES: %.3f',success);
%     disp(d);
%     d = sprintf('(1+1)-ES: %.3f',success1);
%     disp(d);
% 
% end