function [X, FVAL] = Deseason(temperature, seasonFunction, days, guess, options) 
%DESEASON of temperature series
%   Used for the deseasoning of temperature series. Requires input of the
%   temperature, anonymous season function, and fmincon options. Will be
%   extended to correctly handle dynamic time periods.

[X, FVAL] = fmincon(@(x) sum((temperature - ...
    seasonFunction(x, 0:days-1)).^2), guess, ...
    [], [], [], [], [], [], [], options);
end

