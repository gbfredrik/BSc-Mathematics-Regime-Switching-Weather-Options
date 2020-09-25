function [log_f] = log_f_GH(x, lambda, alpha, beta, delta, mu)
%LOG_F_GH Summary of this function goes here
%   Detailed explanation goes here

log_f = lambda/2 * log(alpha^2 - beta^2) ... % TODO: Kontrollera om det blir lambda/2 här pga sqrt(a^2-b^2) !!!
    - (1 / 2 * log(2 * pi) + lambda * log(delta) + (lambda - 1 / 2) * log(alpha) + log(besselk(lambda, delta * sqrt(alpha^2 - beta^2)))) ...
    + (lambda / 2 - 1 / 4) * log(delta^2 + (x - mu)^2) ...
    + log(besselk(lambda - 1/2, alpha * sqrt(delta^2 + (x - mu)^2))) ...
    + beta * (x - mu);


end

