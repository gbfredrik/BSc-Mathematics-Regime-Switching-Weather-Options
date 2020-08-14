function [Set_f, Theta_f, iter_f] = EM(Set, Theta, maxIter, printToggle)

% Iteration settings
iter_f = 0;
delta_Theta = 1;
eps = 1e-6;

% Initial assignment
Set_f = Set;
Theta_f = Theta;

while (delta_Theta > eps)
    iter_f = iter_f + 1;
    
    % Expectation step
    Set_f = EM_Expectation();
    
    % Maximization step
    Theta_f = EM_Maximization();
    
    deltaTheta = norm(Theta_f - Theta, 2);
    Theta = Theta_f; % Remember old optima
    
    if printToggle
        fprintf("Iteration: %d, deltaTheta: %d",...
            iter_f, deltaTheta);
    end
    
    if (iter_f >= maxIter)
        break
    end
end


end