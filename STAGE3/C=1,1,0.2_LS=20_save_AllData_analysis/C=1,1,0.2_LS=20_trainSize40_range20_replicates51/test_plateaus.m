% Print objective function value

% s = 3:5;  % strategy index
% r = ;   % replicate index 

% Plot all plateaus runs 

for p = 11:1:21
    f_num = 800;        % the # of obejctive functions want to return value
    % index of runs plateaus
    r_array3 = find(squeeze(T_f7(3,p,:) == 9999))'; 
    r_array4 = find(squeeze(T_f7(4,p,:) == 9999))';
    r_array5 = find(squeeze(T_f7(5,p,:) == 9999))';
    % objective function value of runs plateaus with objFunCall = f_num
    fx_array3 = squeeze(f_x_f7(3,p,r_array3,f_num))'
    fx_array4 = squeeze(f_x_f7(4,p,r_array4,f_num))'
    fx_array5 = squeeze(f_x_f7(5,p,r_array5,f_num))'


    % r_array = r_array3;
    for j = 1:3
        if j==1 
            r_array = r_array3; 
            s = 3;
        elseif j == 2
            r_array = r_array4; 
            s = 4;
        elseif j==3
            r_array = r_array5; 
            s = 5;
        end
        if length(r_array) ~= 0
            for i = 1:length(r_array)
                f = 7;
                r = r_array(i);
                f_x_array = squeeze(f_x_f7(s,p,r,:));
                sigma_array = squeeze(sigma_f7(s,p,r,:));
                sigma_star_array = squeeze(sigma_star_f7(s,p,r,:));
                success_rate = squeeze(success_f7(s,p,r));
                eval_rate = squeeze(eval_rate_f7(s,p,r));
                four_probs = squeeze(four_prob_f7(s,p,r,:));

                T_range = 1:9999;

                figure(p)

                subplot(1,3,1)
                hold on;
                plot(T_range,f_x_array(T_range));hold on;
                xlabel('function calls','FontSize',15);
                ylabel('function value','FontSize',15);
                set(gca, 'YScale', 'log');

                subplot(1,3,2)
                hold on;
                plot(T_range,sigma_array(T_range));hold on;
                xlabel('function calls','FontSize',15);
                ylabel('step size','FontSize',15);
                set(gca, 'YScale', 'log');

                subplot(1,3,3)
                hold on;
                plot(T_range,sigma_star_array(T_range));hold on;
                xlabel('function calls','FontSize',15);
                ylabel('normalized step size','FontSize',15);
                set(gca, 'YScale', 'log');
            end
        end
    end


    end


