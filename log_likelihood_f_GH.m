function [likelihood] = log_likelihood_f_GH(x, lambda, alpha, beta, delta, mu)
%LOG_SUM_F_GH Summary of this function goes here
%   Detailed explanation goes here

likelihood = 0;
%data = [];
for x_i = x'
    %data = [ data f_GH(x_i, lambda, alpha, beta, delta, mu)];
    likelihood = likelihood + log_f_GH(x_i, lambda, alpha, beta, delta, mu);
end

%[lambda, alpha, beta, delta, mu]
end

