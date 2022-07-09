function [phi_i] = generate_phi_i(B_hat_intercept_removed,index_phi)
%GENERATE_PHI_I Summary of this function goes here
%   Detailed explanation goes here
%{
B_hat_intercept_removed - Kxp matrix - contains only the matrix coefficients
index_phi - scalar - index of phi that needs to be calculated, essentially this represents h in phi_h

B_hat_intercept_removed = [.5,.1,0,0;.4,.5,.25,0]

%}


K = size(B_hat_intercept_removed,1);
p = size(B_hat_intercept_removed,2)/K; % calculating the VAR p, what is the lag order

phi_0 = eye(K);

phi_cummulative = phi_0;


if index_phi==0
    phi_i = phi_0;
else
    for i = 1:index_phi

        sum = 0;
        for j = 1:i

            % Checking if the h>p then the coefficient A = 0
            if j>p
                A = zeros(K,K);
            else
                index_A_i = K*(j-1) + 1 : j*K;
                A = B_hat_intercept_removed(:,index_A_i);
            end


            index_phi_i_minus_j = K*(i-j)+1:K*(i-j+1);
            
            sum = sum + phi_cummulative(:,index_phi_i_minus_j) * A;

        end

        phi_cummulative = [phi_cummulative, sum];



    end



    phi_i = phi_cummulative(:,K*(index_phi)+1:K*(index_phi+1));
end
    
  

 




end

