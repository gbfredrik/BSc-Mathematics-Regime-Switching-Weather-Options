function [X, FVAL] = Deseason(temperature, seasonFunction, options) 
%DESEASON of temperature series
%   Used for the deseasoning of temperature series. Requires input of the
%   temperature, anonymous season function, and fmincon options. Will be
%   extended to correctly handle dynamic time periods.

[X, FVAL] = fmincon(@(x) sum((temperature - ...
    seasonFunction(x(1), x(2), x(3), x(4), 0:3650)).^2), ... % TODO Change time periods to dynamic
    [18;0.005;5;0], [], [], [], [], [], [], [], options);
end

