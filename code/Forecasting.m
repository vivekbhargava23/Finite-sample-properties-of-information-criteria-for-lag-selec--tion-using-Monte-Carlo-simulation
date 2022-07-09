function [Y_forecasted,MSE_y_hat] = Forecasting(Y_t,p,h,intercept)
%FORECASTING Summary of this function goes here

%{
Y_t -  k x (T+p) - matrix of y values
p - scalar - lagged order of VAR p 
h - scalar - no. of lags
intercept - scalar - whether there is an intercept or not


%}

T1 = size(Y_t,2);

[B_hat,residual_covariance_hat, t_ratio,ZZ_prime,ZZ_prime_by_T] = VAR_est(Y_t',p,intercept);

K = size(B_hat,1);

mu = zeros(K,1);

if intercept ==1
    mu = B_hat(:,1);
    B_hat_intercept_removed = B_hat(:,2:end);  % separating the Mu column from the B matrix
    
else 
    B_hat_intercept_removed = B_hat;
    mu = zeros(K,1);
end



% Forecasting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = T1 +1 : T1 + h
    sum = zeros(K,1);
    for j = 1:p
        
        index = K*(j-1) + 1 : j*K ;
        
        A = B_hat_intercept_removed(:,index); % these A are one of the coefficient matrices of VAR(p)
        
        sum = sum + A* Y_t(:,i-j);
        
    end
    
    Y_t(:,i) = mu + sum;
        
 
end

Y_forecasted = Y_t(:,T1 +1:T1 + h);
%Y_t = Y_t(:,1:T1);  % converting Y_t back to its original form by removing the forecasted values 


% Calculating MSE y_h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum = 0;
for i = 0:h-1
    
    phi_i = generate_phi_i(B_hat_intercept_removed,i);
    
    sum = sum + phi_i * residual_covariance_hat * phi_i';
    
end

MSE_y_h = sum;



% Creating B - matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = zeros(K*p+1,K*p+1);
B(1,1) = 1;
B(2:2+K-1,:) = [mu B_hat_intercept_removed];


if p>1
    I_k_p_minus_1 = eye(K*(p-1));
    B(2+K:end,2:K*p + 1 - K) = I_k_p_minus_1;
end


if intercept ==0
    B = B(2:end,2:end)
    
end


% Calculating omega
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega_h = 0;
for i = 0:h-1
    
    phi_i = generate_phi_i(B_hat_intercept_removed,i); % this gives us phi_i
    
    for j = 0:h-1
        
        phi_j = generate_phi_i(B_hat_intercept_removed,j); 
        
        First_inside_term = (B')^(h-1-i);
        Second_inside_term = inv(ZZ_prime_by_T);
        Third_inside_term = (B)^(h-1-j);
        Fourth_inside_term = ZZ_prime_by_T;
        
        trace_value  = trace(First_inside_term * Second_inside_term * Third_inside_term * Fourth_inside_term);
        
        omega_h = omega_h + trace_value * phi_i * residual_covariance_hat * phi_j';
        
        
    end
    
end


% Calculating MSE_y_hat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSE_y_hat = MSE_y_h + omega_h/(T1-p); % (T1-p) to account for pre-sample values




        
   








end

