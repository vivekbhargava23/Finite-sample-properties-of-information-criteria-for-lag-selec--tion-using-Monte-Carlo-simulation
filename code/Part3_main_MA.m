%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Part 3 - MA case with 
lambda1_star = 0.2
lambda2_star = 0.4

& 

lambda1_star = 0.6
lambda2_star = 0.8


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

clear all;
close all;
clc;

rng('default')


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
T_range = [50 100 200];
h_range = [1,2,4];
M = 5000 % 00
p_max = 8;
K = 2;
B = [1,0;0.5,1];
no_of_info_criteria = 4; % [FPE AIC HQ SC]

mu = zeros(K,1);
var_covar = eye(K);
var_covar_ut = B * B';

% values for plots
fontsize_subplots = 12;


% 1.f Repeating for MA(1) process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case 3 - MA(1) with Inital Lambda values
lambda1_star = 0.2
lambda2_star = 0.4

for t = 1:size(T_range,2)
    
    T = T_range(t)
    
    sum_MSPE_h1_MA1 = zeros(1,no_of_info_criteria);
    sum_MSPE_h2_MA1 = zeros(1,no_of_info_criteria);
    sum_MSPE_h4_MA1 = zeros(1,no_of_info_criteria); % initiallizing all the values as 0
    
    
    for i = 1:M
        
        % for each 5000 replications
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Z = normrnd (0,1, [K K]);
        X = Z * Z';
        [U,S,V] = svd(X);
        lambda_star = diag([lambda1_star lambda2_star]);
        
        A = V * lambda_star * V';
        
        
        n = T + p_max + 4 + 1; % adding 1 because now we have a MA1 process
        
        epsilon_t_MA1 = mvnrnd(mu,var_covar,n)';
        
        u_t_MA1 = B * epsilon_t_MA1;
        
        % Generating the series with T+p_max+4 observations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % y_t(:,1) = zeros(K,1);
        
        for ii = 1:n-1   % the series starts from 2 since its a VAR 1 process and thats why we need to add 1 to n so that we have T + pmax + 4
            
            y_t_MA1(:,ii) = u_t_MA1(:,ii+1) + A * u_t_MA1(:,ii); 
            
        end
        
        
        % Generating the series with T+p_max observations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Y_t_MA1 = y_t_MA1(:,1:T + p_max); % seggregating the Y_t into time period till T and T+4
        Y_T_plus_h_MA1 = y_t_MA1(:,T + p_max + 1:end); % seggregating the Y_t into time period till T and T+4
        
        % Calculating the lag order for the computed time series
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        intercept = 1;
        [~,~,Lag_order] = lag_order_v2(Y_t_MA1,p_max,intercept);
        
        Lag_order_frequency_MA1(i,:,t) = Lag_order;
        
        
        % Calculating the h-step ahead forecast
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for pp = 1:size(Lag_order,2)
            
            p = Lag_order(pp); 
            h=4;

            if p == 0

                [B_hat,residual_covariance_hat, t_ratio]= VAR_est(Y_t_MA1',p,intercept); % since p=0 B_hat will be just the mu value
                Y_forecasted = B_hat .* ones(2,4); % for any value of h the forecast will be same when p=0.

            else

                [Y_forecasted,~] = Forecasting(Y_t_MA1,p,h,1); % this gives forecasted value for h= 1 2 3 and 4

            end


            % Calculating the MSE of true DGP
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            MSE_y_h_MA1 = var_covar_ut + A * var_covar_ut * A'; % The MSE matrix in case of MA1 is different than VAR1 process


            % for h=1
            diff_h1 = Y_T_plus_h_MA1(:,1) - Y_forecasted(:,1);
            sum_MSPE_h1_MA1(pp) = sum_MSPE_h1_MA1(pp) + diff_h1' * inv(MSE_y_h_MA1) * diff_h1;

            % for h = 2
            diff_h2 = Y_T_plus_h_MA1(:,2) - Y_forecasted(:,2);
            sum_MSPE_h2_MA1(pp) = sum_MSPE_h2_MA1(pp) + diff_h2' * inv(MSE_y_h_MA1) * diff_h2;

            % for h=4
            diff_h4 = Y_T_plus_h_MA1(:,4) - Y_forecasted(:,4);
            sum_MSPE_h4_MA1(pp) = sum_MSPE_h4_MA1(pp) + diff_h4' * inv(MSE_y_h_MA1) * diff_h4;
        end
        
        
    end
    
    
    MSPE_h_MA1(1,:,t) = sum_MSPE_h1_MA1 / M;
    MSPE_h_MA1(2,:,t) = sum_MSPE_h2_MA1 / M;
    MSPE_h_MA1(3,:,t) = sum_MSPE_h4_MA1 / M;
    
     
    % Calculating frequency
    for l = 0:p_max
        frequency_table_MA1(l+1,:,t) = sum(Lag_order_frequency_MA1(:,:,t)==l);
        
    end
     
    
end

Relative_frequency_MA1 = frequency_table_MA1/M;

% MSPE Table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T3_MSPE_MA1_FPE_set1 = array2table([MSPE_h_MA1(:,1,1) MSPE_h_MA1(:,1,2) MSPE_h_MA1(:,1,3)],'RowNames',{'h=1','h=2','h=4'},'VariableNames',{'T=50','T=100','T=200'})
T3_MSPE_MA1_AIC_set1 = array2table([MSPE_h_MA1(:,2,1) MSPE_h_MA1(:,2,2) MSPE_h_MA1(:,2,3)],'RowNames',{'h=1','h=2','h=4'},'VariableNames',{'T=50','T=100','T=200'})
T3_MSPE_MA1_HQ_set1 = array2table([MSPE_h_MA1(:,3,1) MSPE_h_MA1(:,3,2) MSPE_h_MA1(:,3,3)],'RowNames',{'h=1','h=2','h=4'},'VariableNames',{'T=50','T=100','T=200'})
T3_MSPE_MA1_SC_set1 = array2table([MSPE_h_MA1(:,4,1) MSPE_h_MA1(:,4,2) MSPE_h_MA1(:,4,3)],'RowNames',{'h=1','h=2','h=4'},'VariableNames',{'T=50','T=100','T=200'})


% Relative frequency chart
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
column_names = {'p=0','p=1','p=2','p=3','p=4','p=5','p=6','p=7','p=8',};
row_names_time = {'T=50','T=100','T=200'};

Rel_freq_MA1_FPE = array2table([Relative_frequency_MA1(:,1,1)';Relative_frequency_MA1(:,1,2)';Relative_frequency_MA1(:,1,3)'],'VariableNames',column_names,'RowNames',row_names_time)
Rel_freq_MA1_AIC = array2table([Relative_frequency_MA1(:,2,1)';Relative_frequency_MA1(:,2,2)';Relative_frequency_MA1(:,2,3)'],'VariableNames',column_names,'RowNames',row_names_time)
Rel_freq_MA1_HQ = array2table([Relative_frequency_MA1(:,3,1)';Relative_frequency_MA1(:,3,2)';Relative_frequency_MA1(:,3,3)'],'VariableNames',column_names,'RowNames',row_names_time)
Rel_freq_MA1_SC = array2table([Relative_frequency_MA1(:,4,1)';Relative_frequency_MA1(:,4,2)';Relative_frequency_MA1(:,4,3)'],'VariableNames',column_names,'RowNames',row_names_time)


%{
% FPE
figure
subplot(1,3,1)
hist(Lag_order_frequency_MA1(:,1,1))
ylim([0 M])
%xlim([0 5])
title('T=50')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,2)
hist(Lag_order_frequency_MA1(:,1,2))
ylim([0 M])
%xlim([0 5])
title('T=100')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,3)
hist(Lag_order_frequency_MA1(:,1,3))
ylim([0 M])
%xlim([0 5])
title('T=200')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

sgtitle(['\fontsize{18}Case 3 - MA(1) - Lag Orders FPE ', newline, '\fontsize{15}(\lambda^*_1 = 0.2  \lambda^*_2 = 0.4)']  ,'fontsize' , 18,'FontName' , 'Times')


% AIC
figure
subplot(1,3,1)
hist(Lag_order_frequency_MA1(:,1,1))
ylim([0 M])
%xlim([0 5])
title('T=50')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,2)
hist(Lag_order_frequency_MA1(:,1,2))
ylim([0 M])
%xlim([0 5])
title('T=100')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,3)
hist(Lag_order_frequency_MA1(:,1,3))
ylim([0 M])
%xlim([0 5])
title('T=200')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

sgtitle(['\fontsize{18}Case 3 - MA1(1) - Lag Orders AIC ', newline, '\fontsize{15}(\lambda^*_1 = 0.2  \lambda^*_2 = 0.4)']  ,'fontsize' , 18,'FontName' , 'Times')





% HQ
figure
subplot(1,3,1)
hist(Lag_order_frequency_MA1(:,3,1))
ylim([0 M])
%xlim([0 5])
title('T=50')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,2)
hist(Lag_order_frequency_MA1(:,3,2))
ylim([0 M])
%xlim([0 5])
title('T=100')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,3)
hist(Lag_order_frequency_MA1(:,3,3))
ylim([0 M])
%xlim([0 5])
title('T=200')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

sgtitle(['\fontsize{18}Case 3 - MA(1) - Lag Orders HQ ', newline, '\fontsize{15}(\lambda^*_1 = 0.2  \lambda^*_2 = 0.4)']  ,'fontsize' , 18,'FontName' , 'Times')



% SC
figure
subplot(1,3,1)
hist(Lag_order_frequency_MA1(:,4,1),8)
ylim([0 M])
%xlim([0 5])
title('T=50')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,2)
hist(Lag_order_frequency_MA1(:,4,2),8)
ylim([0 M])
%xlim([0 5])
title('T=100')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,3)
hist(Lag_order_frequency_MA1(:,4,3),8)
ylim([0 M])
%xlim([0 5])
title('T=200')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

sgtitle(['\fontsize{18}Case 3 - MA(1) - Lag Orders SC ', newline, '\fontsize{15}(\lambda^*_1 = 0.2  \lambda^*_2 = 0.4)']  ,'fontsize' , 18,'FontName' , 'Times')
%}



% 1.f Repeating for MA(1) process with new Lambda values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Case 4 - MA(1) with New Lambda values
lambda1_star = 0.6
lambda2_star = 0.8

%

for t = 1:size(T_range,2)
    
    T = T_range(t)
    
    sum_MSPE_h1_MA1 = zeros(1,no_of_info_criteria);
    sum_MSPE_h2_MA1 = zeros(1,no_of_info_criteria);
    sum_MSPE_h4_MA1 = zeros(1,no_of_info_criteria); % initiallizing all the values as 0
    
    
    for i = 1:M
        
        % for each 5000 replications
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Z = normrnd (0,1, [K K]);
        X = Z * Z';
        [U,S,V] = svd(X);
        lambda_star = diag([lambda1_star lambda2_star]);
        
        A = V * lambda_star * V';
        
        n = T + p_max + 4 + 1; % adding 1 because now we have a MA1 process
        
        epsilon_t_MA1 = mvnrnd(mu,var_covar,n)';
        
        u_t_MA1 = B * epsilon_t_MA1;
        
        % Generating the series with T+p_max+4 observations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % y_t(:,1) = zeros(K,1);
        
        for ii = 1:n-1   % the series starts from 2 since its a VAR 1 process and thats why we need to add 1 to n so that we have T + pmax + 4
            
            y_t_MA1(:,ii) = u_t_MA1(:,ii+1) + A * u_t_MA1(:,ii); 
            
        end
        
        
        % Generating the series with T+p_max observations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Y_t_MA1 = y_t_MA1(:,1:T + p_max); % seggregating the Y_t into time period till T and T+4
        Y_T_plus_h_MA1 = y_t_MA1(:,T + p_max + 1:end); % seggregating the Y_t into time period till T and T+4
        
        
        
        
        % Calculating the lag order for the computed time series
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        intercept = 1;
        [~,~,Lag_order] = lag_order_v2(Y_t_MA1,p_max,intercept);
        
        Lag_order_frequency_MA1_second(i,:,t) = Lag_order;
        
        
        % Calculating the h-step ahead forecast
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for pp = 1:size(Lag_order,2)
            
            p = Lag_order(pp); 
            h=4;

            if p == 0

                [B_hat,residual_covariance_hat, t_ratio]= VAR_est(Y_t_MA1',p,intercept); % since p=0 B_hat will be just the mu value
                Y_forecasted = B_hat .* ones(2,4); % for any value of h the forecast will be same when p=0.

            else

                [Y_forecasted,~] = Forecasting(Y_t_MA1,p,h,1); % this gives forecasted value for h= 1 2 3 and 4

            end


            % Calculating the MSE of true DGP
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            MSE_y_h_MA1 = var_covar_ut + A * var_covar_ut * A'; % The MSE matrix in case of MA1 is different than VAR1 process


            % for h=1
            diff_h1 = Y_T_plus_h_MA1(:,1) - Y_forecasted(:,1);
            sum_MSPE_h1_MA1(pp) = sum_MSPE_h1_MA1(pp) + diff_h1' * inv(MSE_y_h_MA1) * diff_h1;

            % for h = 2
            diff_h2 = Y_T_plus_h_MA1(:,2) - Y_forecasted(:,2);
            sum_MSPE_h2_MA1(pp) = sum_MSPE_h2_MA1(pp) + diff_h2' * inv(MSE_y_h_MA1) * diff_h2;

            % for h=4
            diff_h4 = Y_T_plus_h_MA1(:,4) - Y_forecasted(:,4);
            sum_MSPE_h4_MA1(pp) = sum_MSPE_h4_MA1(pp) + diff_h4' * inv(MSE_y_h_MA1) * diff_h4;
        end
        
        
    end
    MSPE_h_MA1_new_lambda_star(1,:,t) = sum_MSPE_h1_MA1 / M;
    MSPE_h_MA1_new_lambda_star(2,:,t) = sum_MSPE_h2_MA1 / M;
    MSPE_h_MA1_new_lambda_star(3,:,t) = sum_MSPE_h4_MA1 / M;
    

    for l = 0:p_max
        frequency_table_MA1_new_lambda(l+1,:,t) = sum(Lag_order_frequency_MA1_second(:,:,t)==l);
        
    end
     
    
    
end

Relative_frequency_MA1_new_lambda = frequency_table_MA1_new_lambda/M;

% MSPE Table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T4_MSPE_MA1_FPE_set2 = array2table([MSPE_h_MA1_new_lambda_star(:,1,1) MSPE_h_MA1_new_lambda_star(:,1,2) MSPE_h_MA1_new_lambda_star(:,1,3)],'RowNames',{'h=1','h=2','h=4'},'VariableNames',{'T=50','T=100','T=200'})
T4_MSPE_MA1_AIC_set2 = array2table([MSPE_h_MA1_new_lambda_star(:,2,1) MSPE_h_MA1_new_lambda_star(:,2,2) MSPE_h_MA1_new_lambda_star(:,2,3)],'RowNames',{'h=1','h=2','h=4'},'VariableNames',{'T=50','T=100','T=200'})
T4_MSPE_MA1_HQ_set2 = array2table([MSPE_h_MA1_new_lambda_star(:,3,1) MSPE_h_MA1_new_lambda_star(:,3,2) MSPE_h_MA1_new_lambda_star(:,3,3)],'RowNames',{'h=1','h=2','h=4'},'VariableNames',{'T=50','T=100','T=200'})
T4_MSPE_MA1_SC_set2 = array2table([MSPE_h_MA1_new_lambda_star(:,4,1) MSPE_h_MA1_new_lambda_star(:,4,2) MSPE_h_MA1_new_lambda_star(:,4,3)],'RowNames',{'h=1','h=2','h=4'},'VariableNames',{'T=50','T=100','T=200'})


% Relative frequency chart
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
column_names = {'p=0','p=1','p=2','p=3','p=4','p=5','p=6','p=7','p=8',};
row_names_time = {'T=50','T=100','T=200'};

Rel_freq_MA1_FPE_new_lambda = array2table([Relative_frequency_MA1_new_lambda(:,1,1)';Relative_frequency_MA1_new_lambda(:,1,2)';Relative_frequency_MA1_new_lambda(:,1,3)'],'VariableNames',column_names,'RowNames',row_names_time)
Rel_freq_MA1_AIC_new_lambda = array2table([Relative_frequency_MA1_new_lambda(:,2,1)';Relative_frequency_MA1_new_lambda(:,2,2)';Relative_frequency_MA1_new_lambda(:,2,3)'],'VariableNames',column_names,'RowNames',row_names_time)
Rel_freq_MA1_HQ_new_lambda = array2table([Relative_frequency_MA1_new_lambda(:,3,1)';Relative_frequency_MA1_new_lambda(:,3,2)';Relative_frequency_MA1_new_lambda(:,3,3)'],'VariableNames',column_names,'RowNames',row_names_time)
Rel_freq_MA1_SC_new_lambda = array2table([Relative_frequency_MA1_new_lambda(:,4,1)';Relative_frequency_MA1_new_lambda(:,4,2)';Relative_frequency_MA1_new_lambda(:,4,3)'],'VariableNames',column_names,'RowNames',row_names_time)



%{
% FPE
figure
subplot(1,3,1)
hist(Lag_order_frequency_MA1_second(:,1,1))
ylim([0 M])
%xlim([0 5])
title('T=50')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,2)
hist(Lag_order_frequency_MA1_second(:,1,2))
ylim([0 M])
%xlim([0 5])
title('T=100')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,3)
hist(Lag_order_frequency_MA1_second(:,1,3))
ylim([0 M])
%xlim([0 5])
title('T=200')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

sgtitle(['\fontsize{18}Case 4 - MA(1) - Lag Orders FPE ', newline, '\fontsize{15}(\lambda^*_1 = 0.6  \lambda^*_2 = 0.8)']  ,'fontsize' , 18,'FontName' , 'Times')


% AIC
figure
subplot(1,3,1)
hist(Lag_order_frequency_MA1_second(:,1,1))
ylim([0 M])
%xlim([0 5])
title('T=50')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,2)
hist(Lag_order_frequency_MA1_second(:,1,2))
ylim([0 M])
%xlim([0 5])
title('T=100')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,3)
hist(Lag_order_frequency_MA1_second(:,1,3))
ylim([0 M])
%xlim([0 5])
title('T=200')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

sgtitle(['\fontsize{18}Case 4 - MA1(1) - Lag Orders AIC ', newline, '\fontsize{15}(\lambda^*_1 = 0.6  \lambda^*_2 = 0.8)']  ,'fontsize' , 18,'FontName' , 'Times')





% HQ
figure
subplot(1,3,1)
hist(Lag_order_frequency_MA1_second(:,3,1))
ylim([0 M])
%xlim([0 5])
title('T=50')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,2)
hist(Lag_order_frequency_MA1_second(:,3,2))
ylim([0 M])
%xlim([0 5])
title('T=100')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,3)
hist(Lag_order_frequency_MA1_second(:,3,3))
ylim([0 M])
%xlim([0 5])
title('T=200')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

sgtitle(['\fontsize{18}Case 4 - MA(1) - Lag Orders HQ ', newline, '\fontsize{15}(\lambda^*_1 = 0.6  \lambda^*_2 = 0.8)']  ,'fontsize' , 18,'FontName' , 'Times')



% SC
figure
subplot(1,3,1)
hist(Lag_order_frequency_MA1_second(:,4,1),8)
ylim([0 M])
%xlim([0 5])
title('T=50')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,2)
hist(Lag_order_frequency_MA1_second(:,4,2),8)
ylim([0 M])
%xlim([0 5])
title('T=100')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

subplot(1,3,3)
hist(Lag_order_frequency_MA1_second(:,4,3),8)
ylim([0 M])
%xlim([0 5])
title('T=200')
set (gca, 'FontName' , 'Times' , 'fontsize' , fontsize_subplots)

sgtitle(['\fontsize{18}Case 4 - MA(1) - Lag Orders SC ', newline, '\fontsize{15}(\lambda^*_1 = 0.6  \lambda^*_2 = 0.8)']  ,'fontsize' , 18,'FontName' , 'Times')





    

    
    
%}
%}
%}
%}
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
            
            
            
            
            
            
        
    
    
    
    





































