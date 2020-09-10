function [likelihood] = log_likelihood_f_VG(x, lambda, alpha, beta, mu)
%LOG_SUM_F_VG Summary of this function goes here
%   Detailed explanation goes here

likelihood = 0;
%data = [];
for x_i = x'
    %data = [ data f_GH(x_i, lambda, alpha, beta, delta, mu)];
    likelihood = likelihood + log_f_VG(x_i, lambda, alpha, beta, mu);
end

%[lambda, alpha, beta, delta, mu]
end

