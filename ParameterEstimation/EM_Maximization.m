function [Theta, p, Q] = EM_Maximization(T_deseason, Pr1, Pr2, Pr1T, Pr2T, p, theta_base, theta_shift, model)
%EM_MAXIMIZATION Summary of this function goes here
%   Detailed explanation goes here

if (model == 1)
    kappa_n1 = sum(Pr1T(3:end)' .* (T_deseason(2:end-1) ...
        .* (T_deseason(3:end) - T_deseason(2:end-1))) ./ (T_deseason(2:end-1).^2)) ...
        / (sum(Pr1T(3:end)));

    sigma1_n1 = sqrt((sum(Pr1T(3:end)' .* T_deseason(2:end-1).^(-2) ...
        .* (T_deseason(3:end) - (1 + theta_base(1)) .* T_deseason(2:end-1)).^2))...
        / (sum(Pr1T(3:end))));

    mu2_n1 = (sum(Pr2T(3:end)' .* (T_deseason(3:end) - T_deseason(2:end-1)))) ...
        / (sum(Pr2T(3:end)));

    sigma2_n1 = sqrt((sum(Pr2T(3:end)' .* (T_deseason(3:end) - (T_deseason(2:end-1) + theta_shift(1))).^2)) ...
        / (sum(Pr2T(3:end))));
elseif (model == 2)
    
end

% Update transition probabilities
% First, compute forecast probabilities:
Pr1_tminus1 = Pr1(1:end-1) * p(1, 1) + Pr2(1:end-1) * p(2, 1);
Pr2_tminus1 = Pr1(1:end-1) * p(1, 2) + Pr2(1:end-1) * p(2, 2);
% Second, compute the transitions:
p(1, 1) =  sum(Pr1T(2:end) .* p(1, 1) .* ( Pr1(1:end-1) ./ Pr1_tminus1 )) ...
    / sum(Pr1T(1:end-1));
%p(1, 2) = 1 - p(1, 1);
p(1, 2) =  sum(Pr2T(2:end) .* p(1, 2) .* ( Pr1(1:end-1) ./ Pr2_tminus1 )) ...
    / sum(Pr1T(1:end-1));

p(2, 1) =  sum(Pr1T(2:end) .* p(2, 1) .* ( Pr2(1:end-1) ./ Pr1_tminus1 )) ...
    / sum(Pr2T(1:end-1));
p(2, 2) = 1 - p(2, 1);

%p(2, 2) =  sum(Pr2T(2:end) .* p(2, 2) .* ( Pr2(1:end-1) ./ Pr2_tminus1)) ...
%    / sum(Pr2T(1:end-1));

fprintf("p(1, 1): %.3f, p(1, 2): %.3f. p(2, 1): %.3f, p(2, 2): %.3f.\n\n",...
    p(1, 1), p(1, 2), p(2, 1), p(2, 2));

% Todo: Scrap
p1_n1 = 0.98;
p2_n1 = 0.95;

if (model == 1)
    Theta = [kappa_n1, sigma1_n1, p1_n1, mu2_n1, sigma2_n1, p2_n1];
elseif (model == 2)
    
end

Q = [1];
end

