function [series, states] = SimulateHMM(initialProbs, transitionProbs, paramBase, paramShift, model, T, S, paramGH, distGH)
%SIMULATEHMM Simulate a two-state HMM
%   initialProbs: [p_base(t=0) p_shifted(t=0)]
%   transitionProbs: [p_11 p_12; p_21 p_22]
%   paramBase: parameters of base distribution
%   paramShift: parameters of shifted distribution
%   model: 1 = two Wiener processes, 2 = Gyamerah, 3 = Evarest
%   T: number of samples simulated per sequence
%   S: number of sequences simulated


if (length(initialProbs) ~= 2)
    error('initialProbs must take two values');
end

if (size(transitionProbs, 1) ~= 2 && size(transitionProbs, 2) ~= 2)
    error('transitionProbs must be a 2x2 matrix')
end

if (model == 2 && nargin ~= 9)
    error('Model 2 requires paramGH and pdfL to be given')
end

% Initialize output vectors to TxS zeros
states = zeros(T, S);
series = zeros(T, S);

rng('default')

% Loop S simulated sequences
for s = 1:S
    %% Randomize starting state by checking against initialProbs
    start = rand;
    states(1, s) = start <= initialProbs(1);

    if (model == 1)
        if (states(1, s) == 0)
            series(1, s) = paramBase(1) + paramBase(2) * normrnd(0, 1);
        else
            series(1, s) = paramShift(1) + paramShift(2) * normrnd(0, 1);
        end
    end
    % Loop T simulations from given model
    for t = 2:T
        if (states(t-1, s) == 0)
            states(t, s) = 1 - (rand <= transitionProbs(1, 1));
        else
            states(t, s) = rand <= transitionProbs(2, 2);
        end
        
        % Model 1: paramBase(mu, sigma), paramShift(mu, sigma)
        % Model 2: 
        % paramBase(kappa, sigma_1, p_1), paramShift(mu_2, sigma_2, p_2)
        % Model 3: 
        % paramBase(beta, mu_1, sigma_1, p_1), paramShift(mu_2, sigma_2, p_2)
        if (model == 1)
            if (states(t, s) == 0)
                series(t, s) = series(t-1, s) + paramBase(1) + paramBase(2) * normrnd(0, 1);
            else
                series(t, s) = series(t-1, s) + paramShift(1) + paramShift(2) * normrnd(0, 1);
            end
        elseif (model == 2)
            if (states(t, s) == 0)
                series(t, s) = (1 + paramBase(1)) * series(t-1, s) + paramBase(2) * abs(series(t-1, s)) * normrnd(0, 1);
            else
                r = rand;
                f_GH = @(x, lambda, alpha, beta, delta, mu) ...
                    sqrt(alpha^2 - beta^2)^(lambda) ... % TODO: SKA DET VARA SQRT HÄR? Se också i ProjectWrapper.m
                    / (sqrt(2*pi) * delta^lambda * alpha^(lambda - 1/2) * besselk(lambda, delta * sqrt(alpha^2 - beta^2))) ...
                    * (sqrt(delta^2 + (x - mu)^2))^(lambda - 1/2) ...
                    * besselk(lambda - 1/2, alpha * sqrt(delta^2 + (x - mu)^2)) * exp(beta * (x - mu));
                f_VG = @(x, lambda, alpha, beta, mu) ...
                    sqrt(alpha^2 - beta^2)^(2 * lambda) ...
                    * abs(x - mu)^(lambda - 1/2) * besselk(lambda - 1/2, alpha * abs(x - mu)) ...
                    / (sqrt(pi) * gamma(lambda) * (2 * alpha)^(lambda - 1/2)) ...
                    * exp(beta * (x - mu));
                
                %[FVAL, res] = fminsearch(@(x) r - integral(@(y) f_GH(y, paramGH), -inf, x));
                %[X, FVAL, exitflag] = fsolve(@(x) r - integral(@(y) normpdf(y, paramGH), -50, x, 'ArrayValued', true), -20);
                if (distGH == 1)
                    [X, FVAL, exitflag] = fsolve(@(x) r - integral(@(y) f_GH(y, paramGH(1), paramGH(2), paramGH(3), paramGH(4), paramGH(5)), -50, x, 'ArrayValued', true), 0);
                elseif (distGH == 2)
                    [X, FVAL, exitflag] = fsolve(@(x) r - integral(@(y) f_VG(y, paramGH(1), paramGH(2), paramGH(3), paramGH(4)), -50, x, 'ArrayValued', true), 0);
                elseif (distGH == 3)
                    "error"
                    break;
                end
                levyResidual = X;
                test(t,s) = X;
                series(t, s) = series(t-1, s) + paramShift(1) + paramShift(2) * levyResidual;
            end
        elseif (model == 3)
            if (states(t, s) == 0)
                series(t, s) = (1 - paramBase(1)) * series(t-1, s) + paramBase(2) + paramBase(3) * abs(series(t-1, s)) * normrnd(0, 1);
                %drift(t,s) = (1 - paramBase(1)) * series(t-1, s) + paramBase(2);
                %diffus(t,s) = paramBase(3) * abs(series(t-1, s)) * normrnd(0, 1);
                %series(t,s) = drift(t,s) + diffus(t,s);
            else
                series(t, s) = series(t-1, s) + paramShift(1) + paramShift(2) * normrnd(0, 1);
                %drift(t,s) = series(t-1, s) + paramShift(1);
                %diffus(t,s) = paramShift(2) * normrnd(0, 1);
                %series(t,s) = drift(t,s) + diffus(t,s);
            end
            if series(t,s) > 1000
                "hej"
            end
        end
    end
end

figure
subplot(2, 1, 1)
plot(series)
subplot(2, 1, 2)
plot(states)


end
% [series, states] = SimulateHMM([0.94 0.06], [0.92 0.08; 0.30 0.70], [0.01 0.10], [0.025 0.25], 1, 250, 1)
% [series, states] = SimulateHMM([0.7895 0.2105], [0.92 0.08; 0.30 0.70], [0.01 0.10], [0.025 0.25], 1, 250, 1)
