function [Pr1, Pr2, Pr1T, Pr2T, Theta_f, p, Q, iter_f] = EM(T_deseason, Theta, p, model, maxIter, printToggle, useKim)

% Iteration settings
iter_f = 0;
delta_Theta = 1;
eps = 1e-10;

% Initial assignment
Theta_f = Theta;
Pr1T = zeros(1, length(T_deseason));
Pr2T = zeros(1, length(T_deseason));
if (model == 1)
    Pr1T(1) = Theta(3);
    Pr2T(1) = Theta(6);
elseif (model == 2)
    Pr1T(1) = Theta(4);
    Pr2T(1) = Theta(7);
end

% Model 1 Theta: kappa, sigma_1, p_1, mu_2, sigma_2, p_2
% Model 2 Theta: beta, mu_1, sigma_1, p_1, mu_2, sigma_2, p_2

Q = [];

while (delta_Theta > eps)
    iter_f = iter_f + 1;
    
    %% Expectation step
    if (model == 1)
        [Pr1_n, Pr2_n, Pr1T_n, Pr2T_n, Q_n] = EM_Expectation( ...
            T_deseason, ...
            Pr1T, ...
            Pr2T, ...
            p, ...
            Theta_f(iter_f, 1:3), ...
            Theta_f(iter_f, 4:6), ...
            model, ...
            useKim);
    elseif (model == 2)
        [Pr1_n, Pr2_n, Pr1T_n, Pr2T_n, Q_n] = EM_Expectation( ...
            T_deseason, ...
            Pr1T, ...
            Pr2T, ...
            p, ...
            Theta_f(iter_f, 1:4), ...
            Theta_f(iter_f, 5:7), ...
            model, ...
            useKim);
    end

    
    %% Maximization step
    if (model == 1)
        [Theta_f(end+1, :), p_n, Q(end+1, :)] = EM_Maximization( ...
            T_deseason, ...
            Pr1_n, ... % > 0.5?
            Pr2_n, ... % > 0.5?
            Pr1T_n, ... % > 0.5?
            Pr2T_n, ... % > 0.5?
            p, ...
            Theta_f(iter_f, 1:3), ...
            Theta_f(iter_f, 4:6), ...
            model, ...
            useKim);
        
    elseif (model == 2)
        [Theta_f(end+1, :), p_n, Q(end+1, :)] = EM_Maximization( ...
            T_deseason, ...
            Pr1_n, ...
            Pr2_n, ...
            Pr1T_n, ...
            Pr2T_n, ...
            p, ...
            Theta_f(iter_f, 1:4), ...
            Theta_f(iter_f, 5:7), ...
            model, ...
            useKim);
        
    end
    
    deltaTheta = norm(Theta_f(end,:) - Theta, 2);
    Theta = Theta_f(end,:); % Remember old optima
    
    if printToggle
        fprintf("Iteration: %d, deltaTheta: %d, Q: %.3f.\n\n",...
            iter_f, deltaTheta, Q(end, 1));
        fprintf(" ... with parameters:\n");
        disp(Theta_f(end,:));
    end
    
    if (iter_f > 1)
%         if Q(end, 1) < Q(end-1, 1)
%             Theta_f = Theta_f(1:end-1, :);
%             Q = Q(1:end-1, 1);
%             iter_f = iter_f - 1;
%             
%             break
%         else
            Pr1 = Pr1_n;
            Pr2 = Pr2_n;
            Pr1T = Pr1T_n;
            Pr2T = Pr2T_n;
            p = p_n;
%         end
    else
        Pr1 = Pr1_n;
        Pr2 = Pr2_n;
        Pr1T = Pr1T_n;
        Pr2T = Pr2T_n;
        p = p_n;
    end
    
    if (iter_f >= maxIter || deltaTheta < 10^-5)
        break
    end
end

end