function [Pr1, Pr2, Pr1T, Pr2T, Q] = EM_Expectation(T_d, Pr1T, Pr2T, p, theta_base, theta_shift, model, useKim)
%EMEXPECTATION Summary of this function goes here
%   Detailed explanation goes here

% Set starting values
Pr1 = zeros(1, length(T_d));
Pr2 = zeros(1, length(T_d));
Pr1(1) = Pr1T(1);
Pr2(1) = Pr2T(1);

f1 = zeros(1, length(T_d));
f2 = zeros(1, length(T_d));
f1(1) = 0.0;
f2(1) = 0.0;

%figure(); hold on
%plot(0,0, '.');
for t = 2 : length(T_d) % 
    if (model == 1)
        kappa = theta_base(1);
        sigma_1 = theta_base(2);
        mu_2 = theta_shift(1);
        sigma_2 = theta_shift(2);
        
        p_1 = theta_base(3);
        p_2 = theta_shift(3);
        
        f1(t) = 1 / (sigma_1 * abs(T_d(t-1)) * sqrt(2*pi)) ...
            * exp(-1/2 * ((T_d(t) - (1 + kappa) * T_d(t-1)) ...
            / (sigma_1 * abs(T_d(t-1))))^2);
        
        if (f1(t) > 1)
            pause(0.01);
        end
        
        f2(t) = 1 / (sigma_2 * sqrt(2*pi)) ...
            * exp(-1/2 * ((T_d(t) - (T_d(t-1) + mu_2)) ...
            / (sigma_2))^2);
    elseif (model == 2)
        beta = theta_base(1);
        mu_1 = theta_base(2);
        sigma_1 = theta_base(3);
        mu_2 = theta_shift(1);
        sigma_2 = theta_shift(2);
        
        p_1 = theta_base(4);
        p_2 = theta_shift(3);
        
        f1(t) = 1 / (sigma_1 * abs(T_d(t-1)) * sqrt(2*pi)) ...
            * exp(-1/2 * ((T_d(t) - ((1 - beta) * T_d(t-1) + mu_1)) ...
            / (sigma_1 * T_d(t-1)))^2);
        f2(t) = 1 / (sigma_2 * sqrt(2*pi)) ...
            * exp(-1/2 * ((T_d(t) - (T_d(t-1) + mu_2)) ...
            / (sigma_2))^2);
    end
    
    if (~useKim)
        Pr1(t) = p_1 * f1(t) ...
            / (p_1 * f1(t) + p_2 * f2(t));
        Pr2(t) = 1-Pr1(t);
%         Pr2(t) = p_2 * f2(t) ...
%             / (p_1 * f1(t) + p_2 * f2(t));
    else
        % Compute forecast probabilities:
        Pr1_Ftminus1 = Pr1(t-1) * p(1, 1) + Pr2(t-1) * p(2, 1);
        Pr2_Ftminus1 = Pr1(t-1) * p(1, 2) + Pr2(t-1) * p(2, 2);

        % Compute state membership probability for time t
        Pr1(t) = Pr1_Ftminus1 * f1(t) ...
            / (Pr1_Ftminus1 * f1(t) + Pr2_Ftminus1 * f2(t));
        Pr2(t) = Pr2_Ftminus1 * f2(t) ...
            / (Pr1_Ftminus1 * f1(t) + Pr2_Ftminus1 * f2(t));
    end
    %plot(t,f1,'.');
    %pause(0.001)
end

if (~useKim)
    Pr1T = Pr1;
    Pr2T = Pr2;
else
    % Smooth state probabilities using Kim's smoothing algorithm
    Pr1T(length(T_d)) = Pr1(length(T_d));
    Pr2T(length(T_d)) = Pr2(length(T_d));

    for t = length(T_d)-1 : -1 : 1
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

if (model == 1)
    if (~useKim)
        p(1, 1) = p_1;
        p(1, 2) = p_1;
        p(2, 1) = p_2;
        p(2, 2) = p_2;
    end
    
    Q = sum( ...
        Pr1T(2:end)' .* ( ...
            log(p(1, 1)) ...
            - log(sigma_1 * abs(T_d(1:end-1)) .* sqrt(2 * pi)) ...
            - 1/2 .* (T_d(2:end) - (1 + kappa) .* T_d(1:end-1)) ./ (sigma_1 .* T_d(1:end-1)) ...
        ) ...
        + Pr1T(2:end)' .* ( ...
            log(p(2, 1)) ...
            - log(sigma_1 * abs(T_d(1:end-1)) .* sqrt(2 * pi)) ...
            - 1/2 .* (T_d(2:end) - (1 + kappa) .* T_d(1:end-1)) ./ (sigma_1 .* T_d(1:end-1)) ...
        ) ...
        + Pr2T(2:end)' .* ( ...
            log(p(1, 2)) ...
            - log(sigma_2 * sqrt(2 * pi)) ...
            - 1/2 .* (T_d(2:end) - (T_d(1:end-1) + mu_2)) ./ (sigma_2) ...
        ) ...
        + Pr2T(2:end)' .* ( ...
            log(p(2, 2)) ...
            - log(sigma_2 * sqrt(2 * pi)) ...
            - 1/2 .* (T_d(2:end) - (T_d(1:end-1) + mu_2)) ./ (sigma_2) ...
        ) ...
    );
    Q = [Q];
elseif (model == 2)
    Q = -1;
end

% figure()
% subplot(1,2,1)
% plot(sort(f1))
% subplot(1,2,2)
% plot(f2)
% pause(0.01)

end
