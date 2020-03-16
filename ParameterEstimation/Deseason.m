function [TempVal] = Deseason(a0, a1, a2, a3, t) 
%DESEASON Summary of this function goes here
%   Detailed explanation goes here

SeasonFunc = @(a0, a1, a2, a3, t) a0 + a1 * t + a2 * sin(2 * pi / 365 * (t - a3));% + a1 * (sin(2 * pi / 365 * t)).^0.1;
TempVal = SeasonFunc(a0, a1, a2, a3, t);
end

