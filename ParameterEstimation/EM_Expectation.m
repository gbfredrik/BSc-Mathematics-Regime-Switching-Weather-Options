function [Pr1, Pr2, Pr1T, Pr2T] = EM_Expectation(T_deseason, Pr1T, Pr2T, p, theta_base, theta_shift, model, useKim)
%EMEXPECTATION Summary of this function goes here
%   Detailed explanation goes here

% Set starting values
Pr1 = zeros(1, length(T_deseason));
Pr2 = zeros(1, length(T_deseason));
Pr1(1) = Pr1T(1);
Pr2(1) = Pr2T(1);

%figure(); hold on
%plot(0,0, '.');
for t = 2 : length(T_deseason) % 
    if (model == 1)
        kappa = theta_base(1);
        sigma_1 = theta_base(2);
        mu_2 = theta_shift(1);
        sigma_2 = theta_shift(2);
        
        p_1 = theta_base(3);
        p_2 = theta_shift(3);
        
        f1 = 1 / (sigma_1 * abs(T_deseason(t-1)) * sqrt(2*pi)) ...
            * exp(-1/2 * ((T_deseason(t) - (1 + kappa) * T_deseason(t-1)) ...
            / (sigma_1 * T_deseason(t-1)))^2);
        
        f2 = 1 / (sigma_2 * sqrt(2*pi)) ...
            * exp(-1/2 * ((T_deseason(t) - (T_deseason(t-1) + mu_2)) ...
            / (sigma_2))^2);
    elseif (model == 2)
        beta = theta_base(1);
        mu_1 = theta_base(2);
        sigma_1 = theta_base(3);
        mu_2 = theta_shift(1);
        sigma_2 = theta_shift(2);
        
        p_1 = theta_base(4);
        p_2 = theta_shift(3);
        
        f1 = 1 / (sigma_1 * abs(T_deseason(t-1)) * sqrt(2*pi)) ...
            * exp(-1/2 * ((T_deseason(t) - ((1 - beta) * T_deseason(t-1) + mu_1)) ...
            / (sigma_1 * T_deseason(t-1)))^2);
        f2 = 1 / (sigma_2 * sqrt(2*pi)) ...
            * exp(-1/2 * ((T_deseason(t) - (T_deseason(t-1) + mu_2)) ...
            / (sigma_2))^2);
    end
    
    if (~useKim)
        Pr1(t) = p_1 * f1 ...
            / (p_1 * f1 + p_2 * f2);
        Pr2(t) = p_2 * f2 ...
            / (p_1 * f1 + p_2 * f2);
    else
        % Compute forecast probabilities:
        Pr1_Ftminus1 = Pr1(t-1) * p(1, 1) + Pr2(t-1) * p(2, 1);
        Pr2_Ftminus1 = Pr1(t-1) * p(1, 2) + Pr2(t-1) * p(2, 2);

        % Compute state membership probability for time t
        Pr1(t) = Pr1_Ftminus1 * f1 ...
            / (Pr1_Ftminus1 * f1 + Pr2_Ftminus1 * f2);
        Pr2(t) = Pr2_Ftminus1 * f2 ...
            / (Pr1_Ftminus1 * f1 + Pr2_Ftminus1 * f2);
    end
    %plot(t,f1,'.');
    %pause(0.001)
end

if (~useKim)
    Pr1T = Pr1;
    Pr2T = Pr2;
else
    % Smooth state probabilities using Kim's smoothing algorithm
    Pr1T(length(T_deseason)) = Pr1(length(T_deseason));
    Pr2T(length(T_deseason)) = Pr2(length(T_deseason));

    for t = length(T_deseason)-1 : -1 : 1
        % First compute forecast probabilities:
        %Pr1T_tplus1 = Pr1T(t) * p(1, 1) + Pr2T(t) * p(2, 1);
        %Pr2T_tplus1 = Pr1T(t) * p(1, 2) + Pr2T(t) * p(2, 2);
        % ... and:
        Pr1_tplus1 = Pr1(t) * p(1, 1) + Pr2(t) * p(2, 1);
        Pr2_tplus1 = Pr1(t) * p(1, 2) + Pr2(t) * p(2, 2);

        % Smooth probabilities by filtering over the full filtration F_T
        Pr1T(t) = (Pr1(t) * Pr1T(t+1) * p(1, 1)) / (Pr1_tplus1) ...
            + (Pr1(t) * Pr2T(t+1) * p(1, 2)) / (Pr2_tplus1);

        Pr2T(t) = (Pr2(t) * Pr1T(t+1) * p(2, 1)) / (Pr1_tplus1) ...
            + (Pr2(t) * Pr2T(t+1) * p(2, 2)) / (Pr2_tplus1);

        %fprintf("Pr1_tplus1: %.3f, Pr2_tplus1: %.3f. Pr1T(t): %.3f, Pr2T(t): %.3f.\n\n",...
        %    Pr1_tplus1, Pr2_tplus1, Pr1T(t), Pr2T(t));
    end
end

end
